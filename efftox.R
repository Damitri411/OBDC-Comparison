## =============================================================================
## EffTox: fully model-based Bayesian design for two-agent OBDC selection
## Implements the design as formalised in Sec. 5.1 of Mukherjee, Takeda & Wason
## "Optimal Dose Combination Selection in Oncology Trials Using Bayesian Designs"
## (based on Thall & Cook, 2004, Biometrics; Brock et al., 2017, BMC Med Res Methodol)
##
## NOTE: this is a self-contained, dependency-free re-implementation (own random-walk
## Metropolis-within-Gibbs sampler, no JAGS/Stan needed) built directly from the
## bivariate-logit + association-parameter formulation and the Lp-ball desirability
## function given in the manuscript. For a fully validated, peer-reviewed
## implementation of EffTox itself, use the CRAN package `trialr` (Brock, based on
## Stan) or `EffToxDesign` (Brian J. Smith) -- see the chat for links. This script
## is offered as a transparent, from-first-principles alternative for teaching /
## design exploration, and to let you experiment with the exact desirability
## contour and admissibility rules quoted in Sec. 5.1.
## =============================================================================

## -----------------------------------------------------------------------------
## 1. Design configuration
## -----------------------------------------------------------------------------

efftox_config <- function(
    dose1 = c(1, 2, 3, 4),     # standardized (e.g. log-dose) levels of Agent A
    dose2 = c(1, 2, 3, 4),     # standardized levels of Agent B
    thetaE = 0.50, thetaT = 0.35,   # efficacy floor / toxicity ceiling used in admissibility
    alphaT = 0.90, alphaE = 0.90,   # required posterior confidence for admissibility
    xi_EWOC = 0.95,                 # continuous overdose-control cutoff
    # desirability anchors (neutral, too-toxic, too-inefficacious) and Lp power
    p_power = 1,
    # priors: independent normal priors on (beta0..3, gamma0..3); rho on (-1,1) via Fisher z
    prior_beta_mean  = c(-1, 0.5, 0.5, 0.1),  prior_beta_sd  = c(2, 2, 2, 1),
    prior_gamma_mean = c(-1.5, 0.5, 0.5, 0.1), prior_gamma_sd = c(2, 2, 2, 1),
    prior_rho_mean = 0, prior_rho_sd = 1,
    n_iter = 6000, n_burnin = 2000, thin = 2,
    isTITE = FALSE,      # <- if TRUE, use the time-to-event weighted-likelihood
    #    version (continuous accrual, patients contribute a
    #    partial/pending outcome while still within their
    #    assessment window); if FALSE, the usual complete-
    #    outcome, cohort-based design (unchanged behaviour).
    Tmax = 42,           # TITE only: full assessment window (e.g. days) within
    # which an efficacy/toxicity event would be observed
    accrual_rate = 1     # TITE only: patients per unit time (Poisson accrual)
) {
  list(dose1 = dose1, dose2 = dose2, thetaE = thetaE, thetaT = thetaT,
       alphaT = alphaT, alphaE = alphaE, xi_EWOC = xi_EWOC, p_power = p_power,
       prior_beta_mean = prior_beta_mean, prior_beta_sd = prior_beta_sd,
       prior_gamma_mean = prior_gamma_mean, prior_gamma_sd = prior_gamma_sd,
       prior_rho_mean = prior_rho_mean, prior_rho_sd = prior_rho_sd,
       n_iter = n_iter, n_burnin = n_burnin, thin = thin,
       isTITE = isTITE, Tmax = Tmax, accrual_rate = accrual_rate)
}

## -----------------------------------------------------------------------------
## 2. Joint efficacy-toxicity probability model
##    logit(pE) = b0+b1*x1+b2*x2+b3*x1*x2 ;  logit(pT) = g0+g1*x1+g2*x2+g3*x1*x2
##    joint law coupled via a Gumbel-type association (Thall & Cook, 2004):
##      P(E=e,T=t) = pE^e(1-pE)^(1-e) pT^t(1-pT)^(1-t)
##                    + (-1)^(e+t) pE(1-pE) pT(1-pT) * (exp(psi)-1)/(exp(psi)+1)
## -----------------------------------------------------------------------------

inv_logit <- function(z) 1 / (1 + exp(-z))

marginal_probs <- function(beta, gamma, x1, x2) {
  pE <- inv_logit(beta[1]  + beta[2]  * x1 + beta[3]  * x2 + beta[4]  * x1 * x2)
  pT <- inv_logit(gamma[1] + gamma[2] * x1 + gamma[3] * x2 + gamma[4] * x1 * x2)
  list(pE = pE, pT = pT)
}

## log P(E=e, T=t) for one observation, given association psi (log odds ratio).
## `w` is the TITE weight (1 for a fully-observed/matured outcome, in (0,1) for
## a patient still pending within their assessment window); w=1 for every
## observation reduces this to the usual unweighted log-likelihood.
log_joint_prob <- function(e, t, pE, pT, psi, w = 1) {
  base <- pE^e * (1 - pE)^(1 - e) * pT^t * (1 - pT)^(1 - t)
  assoc <- (-1)^(e + t) * pE * (1 - pE) * pT * (1 - pT) * (exp(psi) - 1) / (exp(psi) + 1)
  p <- base + assoc
  w * log(pmax(p, 1e-10))
}

## -----------------------------------------------------------------------------
## 2b. TITE (time-to-event) weighting for continuous accrual
##    Linear TITE weight (Cheung & Chappell, 2000), applied here to EffTox's
##    bivariate (efficacy, toxicity) outcome: a patient still within the Tmax
##    assessment window, with neither endpoint event yet observed, contributes
##    a DOWN-WEIGHTED "no event (yet)" observation to the likelihood, with
##    weight = (time on study)/Tmax. Once an event is observed on either
##    endpoint, or the assessment window completes without one, that patient's
##    outcome is fully weighted (w=1). This is the same weighted-likelihood
##    construction used by TITE-CRM/TITE-BLRM, applied to EffTox's joint
##    efficacy-toxicity likelihood via log_joint_prob()'s `w` argument above.
## -----------------------------------------------------------------------------

## For one patient, given time-on-study `u` (already capped at Tmax) and their
## TRUE event times for E and T (Inf if that event never occurs), return the
## outcome as currently observable (possibly still pending) and its weight.
tite_patient_status <- function(u, event_time_E, event_time_T, Tmax) {
  e_obs <- as.integer(event_time_E <= u)
  t_obs <- as.integer(event_time_T <= u)
  matured <- (u >= Tmax) || (e_obs == 1) || (t_obs == 1)
  w <- if (matured) 1 else max(u / Tmax, 1e-3)
  list(e = e_obs, t = t_obs, w = w)
}

## -----------------------------------------------------------------------------
## 3. Log-posterior (log-likelihood + independent normal priors) for MCMC
##    Parameter vector: theta = (beta0..3, gamma0..3, psi)
## -----------------------------------------------------------------------------

log_posterior <- function(theta, dat, cfg) {
  beta  <- theta[1:4]
  gamma <- theta[5:8]
  psi   <- theta[9]
  
  mp <- marginal_probs(beta, gamma, dat$x1, dat$x2)
  w <- if (!is.null(dat$w)) dat$w else rep(1, length(dat$x1))
  ll <- sum(log_joint_prob(dat$e, dat$t, mp$pE, mp$pT, psi, w))
  
  lp <- sum(dnorm(beta,  cfg$prior_beta_mean,  cfg$prior_beta_sd,  log = TRUE)) +
    sum(dnorm(gamma, cfg$prior_gamma_mean, cfg$prior_gamma_sd, log = TRUE)) +
    dnorm(psi, cfg$prior_rho_mean, cfg$prior_rho_sd, log = TRUE)
  
  ll + lp
}

## Random-walk Metropolis sampler (self-contained; no external MCMC package needed)
run_mcmc <- function(dat, cfg, init = NULL, step_sd = NULL) {
  np <- 9
  if (is.null(init)) init <- c(cfg$prior_beta_mean, cfg$prior_gamma_mean, cfg$prior_rho_mean)
  if (is.null(step_sd)) step_sd <- rep(0.15, np)
  
  theta <- init
  cur_lp <- log_posterior(theta, dat, cfg)
  n_keep <- floor((cfg$n_iter - cfg$n_burnin) / cfg$thin)
  samples <- matrix(NA_real_, n_keep, np)
  acc <- 0; kept <- 0
  
  for (it in 1:cfg$n_iter) {
    prop <- theta + rnorm(np, 0, step_sd)
    prop_lp <- log_posterior(prop, dat, cfg)
    if (log(runif(1)) < (prop_lp - cur_lp)) {
      theta <- prop; cur_lp <- prop_lp; acc <- acc + 1
    }
    if (it > cfg$n_burnin && (it - cfg$n_burnin) %% cfg$thin == 0) {
      kept <- kept + 1
      samples[kept, ] <- theta
    }
  }
  colnames(samples) <- c(paste0("beta", 0:3), paste0("gamma", 0:3), "psi")
  attr(samples, "accept_rate") <- acc / cfg$n_iter
  samples
}

## -----------------------------------------------------------------------------
## 4. Desirability function (Sec. 5.1 Lp-ball contour) and posterior summaries
##    D(pE,pT) = [1 - ( ((pT - thetaT)_+ /(1-thetaT))^p + ((thetaE-pE)_+/thetaE)^p )^{1/p} ]_+
## -----------------------------------------------------------------------------

desirability <- function(pE, pT, thetaE, thetaT, p) {
  a <- pmax(0, pT - thetaT) / (1 - thetaT)
  b <- pmax(0, thetaE - pE) / thetaE
  d <- 1 - (a^p + b^p)^(1 / p)
  pmax(d, 0)
}

## For every dose combination, compute posterior mean pE, pT, expected desirability,
## and the two admissibility tail probabilities required in Sec. 5.1:
##   Pr(pT <= thetaT | data) >= alphaT ,  Pr(pE >= thetaE | data) >= alphaE
## plus the continuous overdose-control exclusion Pr(pT > thetaT | data) > xi_EWOC.
grid_posterior_summary <- function(samples, cfg) {
  J <- length(cfg$dose1); K <- length(cfg$dose2)
  pE_hat <- pT_hat <- Dbar <- matrix(NA_real_, J, K)
  PrT_ok <- PrE_ok <- PrT_overdose <- matrix(NA_real_, J, K)
  
  for (j in 1:J) for (k in 1:K) {
    x1 <- cfg$dose1[j]; x2 <- cfg$dose2[k]
    pE_draws <- inv_logit(samples[, "beta0"]  + samples[, "beta1"]  * x1 +
                            samples[, "beta2"]  * x2 + samples[, "beta3"]  * x1 * x2)
    pT_draws <- inv_logit(samples[, "gamma0"] + samples[, "gamma1"] * x1 +
                            samples[, "gamma2"] * x2 + samples[, "gamma3"] * x1 * x2)
    D_draws <- desirability(pE_draws, pT_draws, cfg$thetaE, cfg$thetaT, cfg$p_power)
    
    pE_hat[j, k] <- mean(pE_draws)
    pT_hat[j, k] <- mean(pT_draws)
    Dbar[j, k]   <- mean(D_draws)
    PrT_ok[j, k] <- mean(pT_draws <= cfg$thetaT)
    PrE_ok[j, k] <- mean(pE_draws >= cfg$thetaE)
    PrT_overdose[j, k] <- mean(pT_draws > cfg$thetaT)
  }
  list(pE_hat = pE_hat, pT_hat = pT_hat, Dbar = Dbar,
       PrT_ok = PrT_ok, PrE_ok = PrE_ok, PrT_overdose = PrT_overdose)
}

## n_at_current: number of patients treated so far AT THE CURRENT DOSE (not the whole
## grid). Following common EffTox practice (Thall & Cook, 2004; Brock et al., 2017),
## the efficacy-acceptability screen is only applied once a minimum amount of local
## evidence has accrued (min_n_screen); before that, a dose is not excluded for
## efficacy alone (though the toxicity/EWOC checks still always apply). This avoids
## the design declaring premature futility off a single unlucky cohort of 2-3
## patients, exactly the situation flagged in Sec. 2.5 ("reproducibility requires
## clear pre-specification of calibration quantities").
admissible_set_efftox <- function(post, cfg, n_by_cell = NULL, min_n_screen = 6) {
  tox_ok <- (post$PrT_ok >= cfg$alphaT) & (post$PrT_overdose <= cfg$xi_EWOC)
  if (is.null(n_by_cell)) {
    eff_ok <- (post$PrE_ok >= cfg$alphaE)
  } else {
    eff_ok <- (post$PrE_ok >= cfg$alphaE) | (n_by_cell < min_n_screen)
  }
  tox_ok & eff_ok
}

## Multiple cells CAN legitimately tie on the ranking value. Two different
## situations use two different policies:
##  - INTERIM dose-transition steps (deciding where to send the NEXT cohort)
##    MUST resolve to exactly one dose -- a trial can only treat patients at
##    one place at a time -- so those still use select_best_with_tiebreak()
##    below (lowest total standardized dose j+k as the tie-break).
##  - The FINAL recommended OBDC is a REPORTED CONCLUSION, not an operational
##    necessity, so it is NOT tie-broken: select_all_tied() returns every
##    admissible cell achieving the max value, and all of them are reported
##    as (jointly) recommended.
select_best_with_tiebreak <- function(idx, values) {
  best_val <- max(values)
  tied <- which(values == best_val)
  n_tied <- length(tied)
  if (n_tied > 1) {
    dose_sum <- idx[tied, 1] + idx[tied, 2]
    tied <- tied[which.min(dose_sum)]
  }
  list(sel = as.numeric(idx[tied, ]), n_tied = n_tied)
}

## Returns EVERY row of `idx` achieving max(values), as a 2-column matrix
## with columns named "j","k" -- no tie-break applied.
select_all_tied <- function(idx, values) {
  best_val <- max(values)
  tied <- which(values == best_val)
  out <- matrix(idx[tied, ], ncol = 2)
  colnames(out) <- c("j", "k")
  out
}

## -----------------------------------------------------------------------------
## 5. Dose-transition rule: two-dimensional neighbourhood, no-skipping
## -----------------------------------------------------------------------------

et_neighbours <- function(j, k, J, K) {
  cand <- rbind(c(j + 1, k), c(j, k + 1), c(j - 1, k), c(j, k - 1), c(j, k))
  cand[cand[, 1] >= 1 & cand[, 1] <= J & cand[, 2] >= 1 & cand[, 2] <= K, , drop = FALSE]
}

next_dose_efftox <- function(current, samples, cfg, n_by_cell = NULL, min_n_screen = 6) {
  j <- current[1]; k <- current[2]
  post <- grid_posterior_summary(samples, cfg)
  A <- admissible_set_efftox(post, cfg, n_by_cell, min_n_screen)
  
  cand <- et_neighbours(j, k, length(cfg$dose1), length(cfg$dose2))
  is_adm <- apply(cand, 1, function(z) A[z[1], z[2]])
  cand_ok <- cand[is_adm, , drop = FALSE]
  
  if (nrow(cand_ok) == 0) {
    return(list(next_dose = NA, stop_trial = TRUE, post = post, admissible = A))
  }
  Dvals <- apply(cand_ok, 1, function(z) post$Dbar[z[1], z[2]])
  pick <- select_best_with_tiebreak(cand_ok, Dvals)
  best <- pick$sel
  list(next_dose = best, stop_trial = FALSE, post = post, admissible = A)
}

## -----------------------------------------------------------------------------
## 6. Trial simulation
## -----------------------------------------------------------------------------

simulate_efftox <- function(cfg, true_pE, true_pT, cohort_size = 3, n_max = 45,
                            start_dose = c(1, 1), min_n_screen = 6, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (isTRUE(cfg$isTITE)) {
    return(simulate_efftox_tite(cfg, true_pE, true_pT, n_max = n_max,
                                start_dose = start_dose, min_n_screen = min_n_screen))
  }
  
  dat <- list(x1 = numeric(0), x2 = numeric(0), e = integer(0), t = integer(0))
  J <- length(cfg$dose1); K <- length(cfg$dose2)
  n_by_cell <- matrix(0, J, K)
  current <- start_dose
  total_n <- 0
  path <- list()
  
  while (total_n < n_max) {
    j <- current[1]; k <- current[2]
    m <- min(cohort_size, n_max - total_n)
    x1 <- cfg$dose1[j]; x2 <- cfg$dose2[k]
    
    e_new <- rbinom(m, 1, true_pE[j, k])
    t_new <- rbinom(m, 1, true_pT[j, k])
    
    dat$x1 <- c(dat$x1, rep(x1, m)); dat$x2 <- c(dat$x2, rep(x2, m))
    dat$e  <- c(dat$e,  e_new);      dat$t  <- c(dat$t,  t_new)
    n_by_cell[j, k] <- n_by_cell[j, k] + m
    total_n <- total_n + m
    
    samples <- run_mcmc(dat, cfg)
    step <- next_dose_efftox(current, samples, cfg, n_by_cell, min_n_screen)
    path[[length(path) + 1]] <- list(dose = current, n = m, e = sum(e_new), t = sum(t_new))
    
    if (step$stop_trial) { message("Trial stopped: no admissible dose pairs remain."); break }
    current <- step$next_dose
  }
  
  final_samples <- run_mcmc(dat, cfg)
  final_post <- grid_posterior_summary(final_samples, cfg)
  A_final <- admissible_set_efftox(final_post, cfg, n_by_cell, min_n_screen)
  
  if (any(A_final)) {
    idx <- which(A_final, arr.ind = TRUE)
    Dvals <- final_post$Dbar[A_final]
    sel <- select_all_tied(idx, Dvals)   # ALL admissible cells tied on max desirability -- no tie-break
  } else {
    sel <- matrix(c(NA_real_, NA_real_), nrow = 1, dimnames = list(NULL, c("j", "k")))
  }
  
  list(data = dat, posterior = final_post, admissible = A_final,
       recommended_OBDC = sel, n_tied_OBDC = nrow(sel), path = path, n_total = total_n, samples = final_samples)
}

## -----------------------------------------------------------------------------
## 6b. TITE trial simulation (continuous accrual, weighted likelihood)
##    Patients arrive one at a time (Poisson process, rate cfg$accrual_rate,
##    matching the accrual_rate convention already used elsewhere in this
##    project, e.g. uTPI-Comb). Each patient's TRUE efficacy/toxicity event
##    times (if the event occurs at all) are drawn uniformly over [0, Tmax] --
##    the standard simulation convention in the TITE-CRM literature. Whenever
##    a new patient is about to be assigned a dose, every earlier patient's
##    outcome is re-weighted using their CURRENT follow-up time (time on study
##    so far), via tite_patient_status(), and the new patient's dose is chosen
##    from the resulting weighted posterior -- so a patient who hasn't yet had
##    time to show toxicity/efficacy still informs the decision, just with
##    less weight than a fully-matured patient. At trial close, the final
##    recommendation uses each patient's FULLY MATURED (unweighted) outcome,
##    matching how a TITE design reports its final result once every patient's
##    assessment window has actually completed.
## -----------------------------------------------------------------------------

simulate_efftox_tite <- function(cfg, true_pE, true_pT, n_max = 45,
                                 start_dose = c(1, 1), min_n_screen = 6) {
  J <- length(cfg$dose1); K <- length(cfg$dose2)
  Tmax <- cfg$Tmax; rate <- cfg$accrual_rate
  
  arrival   <- cumsum(rexp(n_max, rate = rate))   # patient arrival times
  dose_j    <- integer(n_max); dose_k <- integer(n_max)
  ev_time_E <- ev_time_T <- numeric(n_max)        # TRUE event times (Inf = event never occurs)
  n_by_cell <- matrix(0, J, K)
  current <- start_dose
  path <- list()
  
  for (i in seq_len(n_max)) {
    j <- current[1]; k <- current[2]
    dose_j[i] <- j; dose_k[i] <- k
    n_by_cell[j, k] <- n_by_cell[j, k] + 1
    
    ## Simulate whether/when this patient's efficacy and toxicity events
    ## occur, each timed uniformly within the assessment window if it occurs
    ## at all.
    ev_time_E[i] <- if (rbinom(1, 1, true_pE[j, k]) == 1) runif(1, 0, Tmax) else Inf
    ev_time_T[i] <- if (rbinom(1, 1, true_pT[j, k]) == 1) runif(1, 0, Tmax) else Inf
    
    ## Build the weighted dataset as it would look AT THIS PATIENT'S ARRIVAL,
    ## using every earlier patient's follow-up time so far, then pick the dose
    ## for THIS patient from the resulting weighted posterior.
    now <- arrival[i]
    if (i > 1) {
      u <- pmin(now - arrival[1:(i - 1)], Tmax)
      st <- lapply(seq_len(i - 1), function(m) tite_patient_status(u[m], ev_time_E[m], ev_time_T[m], Tmax))
      dat <- list(
        x1 = cfg$dose1[dose_j[1:(i - 1)]], x2 = cfg$dose2[dose_k[1:(i - 1)]],
        e  = vapply(st, function(s) s$e, numeric(1)),
        t  = vapply(st, function(s) s$t, numeric(1)),
        w  = vapply(st, function(s) s$w, numeric(1))
      )
      samples <- run_mcmc(dat, cfg)
      step <- next_dose_efftox(current, samples, cfg, n_by_cell, min_n_screen)
      path[[length(path) + 1]] <- list(dose = current, weighted_n_so_far = round(sum(dat$w), 1))
      if (!step$stop_trial) current <- step$next_dose
    }
  }
  
  ## Final recommendation: by trial close every patient's assessment window
  ## has completed, so use each patient's FULLY MATURED (unweighted) outcome.
  e_final <- as.integer(ev_time_E <= Tmax)
  t_final <- as.integer(ev_time_T <= Tmax)
  dat_final <- list(x1 = cfg$dose1[dose_j], x2 = cfg$dose2[dose_k],
                    e = e_final, t = t_final, w = rep(1, n_max))
  final_samples <- run_mcmc(dat_final, cfg)
  final_post <- grid_posterior_summary(final_samples, cfg)
  A_final <- admissible_set_efftox(final_post, cfg, n_by_cell, min_n_screen)
  
  if (any(A_final)) {
    idx <- which(A_final, arr.ind = TRUE)
    Dvals <- final_post$Dbar[A_final]
    sel <- select_all_tied(idx, Dvals)   # ALL admissible cells tied on max desirability -- no tie-break
  } else {
    sel <- matrix(c(NA_real_, NA_real_), nrow = 1, dimnames = list(NULL, c("j", "k")))
  }
  
  list(data = dat_final, posterior = final_post, admissible = A_final,
       recommended_OBDC = sel, n_tied_OBDC = nrow(sel), path = path, n_total = n_max, samples = final_samples,
       arrival = arrival, event_time_E = ev_time_E, event_time_T = ev_time_T)
}

## -----------------------------------------------------------------------------
## 7. Example
## -----------------------------------------------------------------------------

if (identical(environment(), globalenv())) {
  # NB: alphaT/alphaE control how much posterior confidence is demanded before a
  # dose is declared admissible (Sec. 5.1). With only a handful of patients per
  # cohort, 0.90 confidence can be hard to reach early -- as in any real EffTox
  # trial, this is a genuine feature (the design refuses to certify a dose pair
  # on thin evidence), but for this illustrative run we use more permissive
  # values (0.70) purely so the toy trial has room to move past the first cohort.
  cfg <- efftox_config(dose1 = c(1, 2, 3, 4), dose2 = c(1, 2, 3, 4),
                       thetaE = 0.35, thetaT = 0.35,
                       alphaT = 0.70, alphaE = 0.70,
                       n_iter = 4000, n_burnin = 1500)
  
  true_pT <- matrix(c(0.05, 0.10, 0.15, 0.25,
                      0.10, 0.18, 0.28, 0.38,
                      0.18, 0.28, 0.40, 0.50,
                      0.28, 0.40, 0.52, 0.65), nrow = 4, byrow = TRUE)
  true_pE <- matrix(c(0.08, 0.18, 0.28, 0.30,
                      0.18, 0.32, 0.44, 0.45,
                      0.28, 0.44, 0.52, 0.52,
                      0.30, 0.45, 0.52, 0.52), nrow = 4, byrow = TRUE)
  
  # NB: with small cohorts, an early unlucky run of toxicities/non-responses can
  # legitimately trigger safety/futility stopping -- that is the design working as
  # intended, not a bug. seed 4 below happens to run a full 24-patient trial.
  res <- simulate_efftox(cfg, true_pE, true_pT, cohort_size = 3, n_max = 24, seed = 4)
  if (anyNA(res$recommended_OBDC)) {
    cat("Recommended OBDC: none -- no admissible dose pair.\n")
  } else if (nrow(res$recommended_OBDC) == 1) {
    cat("Recommended OBDC (Agent A level, Agent B level):", res$recommended_OBDC[1, ], "\n")
  } else {
    cat("Recommended OBDC:", nrow(res$recommended_OBDC), "cells tied on desirability (all reported, no tie-break):\n")
    print(res$recommended_OBDC)
  }
  cat("Total patients used:", res$n_total, "\n")
  cat("Posterior mean desirability by dose pair:\n")
  print(round(res$posterior$Dbar, 3))
  
  # --- TITE version: same scenario, but with continuous accrual and a
  # time-to-event weighted likelihood, so decisions don't have to wait for
  # every earlier patient to fully mature before the next one is assigned.
  cfg_tite <- efftox_config(dose1 = c(1, 2, 3, 4), dose2 = c(1, 2, 3, 4),
                            thetaE = 0.35, thetaT = 0.35,
                            alphaT = 0.70, alphaE = 0.70,
                            n_iter = 4000, n_burnin = 1500,
                            isTITE = TRUE, Tmax = 42, accrual_rate = 1)
  set.seed(4)
  res_tite <- simulate_efftox(cfg_tite, true_pE, true_pT, n_max = 24)
  if (anyNA(res_tite$recommended_OBDC)) {
    cat("\n[TITE] Recommended OBDC: none -- no admissible dose pair.\n")
  } else if (nrow(res_tite$recommended_OBDC) == 1) {
    cat("\n[TITE] Recommended OBDC (Agent A level, Agent B level):", res_tite$recommended_OBDC[1, ], "\n")
  } else {
    cat("\n[TITE] Recommended OBDC:", nrow(res_tite$recommended_OBDC), "cells tied on desirability (all reported, no tie-break):\n")
    print(res_tite$recommended_OBDC)
  }
  cat("[TITE] Total patients used:", res_tite$n_total, "\n")
}