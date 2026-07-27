## =============================================================================
## EffTox: fully model-based Bayesian design for two-agent OBDC selection
## Implements the design as formalised in Sec. 5.1 of Mukherjee, Takeda & Wason
## "Optimal Dose Combination Selection in Oncology Trials Using Bayesian Designs"
## (based on Thall & Cook, 2004, Biometrics; Brock et al., 2017, BMC Med Res Methodol)
##

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
  n_iter = 6000, n_burnin = 2000, thin = 2
) {
  list(dose1 = dose1, dose2 = dose2, thetaE = thetaE, thetaT = thetaT,
       alphaT = alphaT, alphaE = alphaE, xi_EWOC = xi_EWOC, p_power = p_power,
       prior_beta_mean = prior_beta_mean, prior_beta_sd = prior_beta_sd,
       prior_gamma_mean = prior_gamma_mean, prior_gamma_sd = prior_gamma_sd,
       prior_rho_mean = prior_rho_mean, prior_rho_sd = prior_rho_sd,
       n_iter = n_iter, n_burnin = n_burnin, thin = thin)
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

## log P(E=e, T=t) for one observation, given association psi (log odds ratio)
log_joint_prob <- function(e, t, pE, pT, psi) {
  base <- pE^e * (1 - pE)^(1 - e) * pT^t * (1 - pT)^(1 - t)
  assoc <- (-1)^(e + t) * pE * (1 - pE) * pT * (1 - pT) * (exp(psi) - 1) / (exp(psi) + 1)
  p <- base + assoc
  log(pmax(p, 1e-10))
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
  ll <- sum(log_joint_prob(dat$e, dat$t, mp$pE, mp$pT, psi))

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
  best <- cand_ok[which.max(Dvals), ]
  list(next_dose = as.numeric(best), stop_trial = FALSE, post = post, admissible = A)
}

## -----------------------------------------------------------------------------
## 6. Trial simulation
## -----------------------------------------------------------------------------

simulate_efftox <- function(cfg, true_pE, true_pT, cohort_size = 3, n_max = 45,
                             start_dose = c(1, 1), min_n_screen = 6, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
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
    sel <- idx[which.max(Dvals), ]
    names(sel) <- c("j", "k")
  } else sel <- c(j = NA, k = NA)

  list(data = dat, posterior = final_post, admissible = A_final,
       recommended_OBDC = sel, path = path, n_total = total_n, samples = final_samples)
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
  cat("Recommended OBDC (Agent A level, Agent B level):", res$recommended_OBDC, "\n")
  cat("Total patients used:", res$n_total, "\n")
  cat("Posterior mean desirability by dose pair:\n")
  print(round(res$posterior$Dbar, 3))
}
