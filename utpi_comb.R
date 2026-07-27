## =============================================================================
## uTPI-Comb: utility-integrated, time-to-event Bayesian design for two-agent
## OBDC selection under staggered/continuous accrual.
## Implements the design as formalised in Sec. 5.4 of Mukherjee, Takeda & Wason
## "Optimal Dose Combination Selection in Oncology Trials Using Bayesian Designs"
## (based on Liang, Yang & Yuan, 2024, J. Univ. Sci. Technol. China 54(12):1206-1216)
##
## NOTE: original re-implementation from the manuscript's equations; no public
## code release for uTPI-Comb (the combination extension) could be located --
## only its single-agent predecessor uTPI has published code (see chat: GitHub
## ruitaolin/uTPI, haoluns/uTPI). The manuscript's own comparative summary
## (Sec. 7.1) states that uTPI-Comb's *decision rule* itself only needs
## conjugate Beta/Dirichlet posteriors for (pE, pT); the hazard model of Sec. 5.4
## is the mechanism used to keep those posteriors current under partial,
## time-to-event follow-up. We implement exactly that: a fractional/TITE-weighted
## Beta-Binomial update (in the spirit of TITE-CRM/TITE-BOIN weighting, which is
## also how Damitri's own TITE-BOIN12 zeros-trick JAGS work handles this), which
## is far more tractable than fitting the full continuous-time PH hazard model in
## Sec. 5.4 while reproducing the same admissible-set / utility decision rule.
## =============================================================================

## -----------------------------------------------------------------------------
## 1. Design configuration
## -----------------------------------------------------------------------------

utpi_comb_config <- function(
  J = 4, K = 4,
  phiT = 0.30, phiE = 0.20,      # overdose / sub-efficacy reference rates
  wE = 1, wT = 1,                # utility weights: U = wE*pE - wT*pT
  aE = 1, bE = 1, aT = 1, bT = 1,# Beta priors for pE, pT
  pi_over  = 0.90,                # cap: Pr(pT > phiT | data_t) < pi_over  -> admissible
  mu_b_source = c(phiE = 0.20, phiT = 0.30), # used to derive benchmark utility mu_b
  assessment_window_T = 42,       # days, DLT assessment window
  assessment_window_E = 56,       # days, efficacy assessment window
  accrual_rate = 1.2,             # patients per week (Poisson process), for staggered accrual
  cohort_size = 1,                # uTPI-Comb allows one-at-a-time / rolling accrual
  n_max = 45,
  start_dose = c(1, 1)
) {
  mu_b <- wE * mu_b_source["phiE"] - wT * mu_b_source["phiT"]
  list(J = J, K = K, phiT = phiT, phiE = phiE, wE = wE, wT = wT,
       aE = aE, bE = bE, aT = aT, bT = bT, pi_over = pi_over, mu_b = as.numeric(mu_b),
       assessment_window_T = assessment_window_T, assessment_window_E = assessment_window_E,
       accrual_rate = accrual_rate, cohort_size = cohort_size, n_max = n_max,
       start_dose = start_dose)
}

## -----------------------------------------------------------------------------
## 2. TITE / fractional-information weighted Beta-Binomial update
##    For a patient enrolled at calendar time s, currently at calendar time t,
##    with assessment window W and (possibly right-censored) event time u <= W:
##      weight = 1                          if event observed (u < W, event = 1)
##      weight = min(1, (t - s) / W)        if event-free so far and follow-up incomplete
##      weight = 1                          if fully followed with no event (u = W, event = 0)
##    "Effective" successes/failures use these fractional weights (TITE-CRM/TITE-BOIN
##    style), giving a real-time posterior at any interim look (Sec 5.4: "near real-
##    time updates of hat{U}_jk from partial data").
## -----------------------------------------------------------------------------

tite_weight <- function(enroll_time, now, window, event_time = NA, event_observed = FALSE) {
  if (event_observed) return(1)
  elapsed <- now - enroll_time
  if (!is.na(event_time) && elapsed >= window) return(1)  # fully followed, no event
  pmin(1, pmax(0, elapsed / window))
}

## Given patient-level data frames for toxicity and efficacy (one row per patient,
## per dose cell), compute the TITE-weighted effective counts and the resulting
## posterior mean / overdose probability for pT, and posterior mean / benchmark
## comparison for the utility U_jk.
utpi_comb_posterior <- function(dat_T, dat_E, cfg, now) {
  J <- cfg$J; K <- cfg$K
  pT_hat <- pE_hat <- U_hat <- matrix(NA_real_, J, K)
  PrT_over <- matrix(NA_real_, J, K)
  PrE_ok <- matrix(NA_real_, J, K)
  PrU_gt_mub <- matrix(NA_real_, J, K)

  for (j in 1:J) for (k in 1:K) {
    dT <- dat_T[dat_T$dose1 == j & dat_T$dose2 == k, , drop = FALSE]
    dE <- dat_E[dat_E$dose1 == j & dat_E$dose2 == k, , drop = FALSE]

    if (nrow(dT) == 0) {
      wT_eff_succ <- 0; wT_eff_n <- 0
    } else {
      wts <- mapply(tite_weight, dT$enroll_time, now, cfg$assessment_window_T,
                    dT$event_time, dT$event_observed)
      wT_eff_succ <- sum(wts * dT$event_observed)
      wT_eff_n    <- sum(wts)
    }
    if (nrow(dE) == 0) {
      wE_eff_succ <- 0; wE_eff_n <- 0
    } else {
      wts <- mapply(tite_weight, dE$enroll_time, now, cfg$assessment_window_E,
                    dE$event_time, dE$event_observed)
      wE_eff_succ <- sum(wts * dE$event_observed)
      wE_eff_n    <- sum(wts)
    }

    aT_post <- cfg$aT + wT_eff_succ; bT_post <- cfg$bT + max(0, wT_eff_n - wT_eff_succ)
    aE_post <- cfg$aE + wE_eff_succ; bE_post <- cfg$bE + max(0, wE_eff_n - wE_eff_succ)

    pT_hat[j, k] <- aT_post / (aT_post + bT_post)
    pE_hat[j, k] <- aE_post / (aE_post + bE_post)
    PrT_over[j, k] <- 1 - pbeta(cfg$phiT, aT_post, bT_post)
    PrE_ok[j, k]   <- 1 - pbeta(cfg$phiE, aE_post, bE_post)   # Pr(pE >= phiE | data)

    U_hat[j, k] <- cfg$wE * pE_hat[j, k] - cfg$wT * pT_hat[j, k]

    # Monte Carlo posterior for Pr(U_jk > mu_b | data): draw independent pT,pE posteriors
    pT_draws <- rbeta(3000, aT_post, bT_post)
    pE_draws <- rbeta(3000, aE_post, bE_post)
    U_draws <- cfg$wE * pE_draws - cfg$wT * pT_draws
    PrU_gt_mub[j, k] <- mean(U_draws > cfg$mu_b)
  }

  list(pT_hat = pT_hat, pE_hat = pE_hat, U_hat = U_hat,
       PrT_over = PrT_over, PrE_ok = PrE_ok, PrU_gt_mub = PrU_gt_mub)
}

## Admissible set (Sec. 5.4 / 7.1): A = { (j,k) : Pr(pT(j,k) > phiT | data) < pi_over }
utpi_comb_admissible <- function(post, cfg) post$PrT_over < cfg$pi_over

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
## 3. Allocation rule: next patient goes to the admissible cell maximising
##    Pr(U_jk > mu_b | data) (Sec. 5.4: "probability-of-superiority allocation")
## -----------------------------------------------------------------------------

utpi_neighbours <- function(j, k, J, K) {
  cand <- rbind(c(j + 1, k), c(j, k + 1), c(j - 1, k), c(j, k - 1), c(j, k))
  cand[cand[, 1] >= 1 & cand[, 1] <= J & cand[, 2] >= 1 & cand[, 2] <= K, , drop = FALSE]
}

utpi_comb_next_dose <- function(current, post, cfg) {
  j <- current[1]; k <- current[2]
  A <- utpi_comb_admissible(post, cfg)
  cand <- utpi_neighbours(j, k, cfg$J, cfg$K)
  is_adm <- apply(cand, 1, function(z) A[z[1], z[2]])
  cand_ok <- cand[is_adm, , drop = FALSE]
  if (nrow(cand_ok) == 0) return(list(next_dose = NA, stop_trial = TRUE))
  score <- apply(cand_ok, 1, function(z) post$PrU_gt_mub[z[1], z[2]])
  pick <- select_best_with_tiebreak(cand_ok, score)
  best <- pick$sel
  list(next_dose = best, stop_trial = FALSE)
}

## -----------------------------------------------------------------------------
## 4. Continuous-accrual trial simulation with staggered enrolment
##    (patients enrol as a Poisson process; dose decisions are made "as of now",
##    using whatever fraction of follow-up each enrolled patient has accrued)
## -----------------------------------------------------------------------------

simulate_utpi_comb <- function(cfg, true_pE, true_pT, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dat_T <- data.frame(dose1 = integer(0), dose2 = integer(0), enroll_time = numeric(0),
                       event_time = numeric(0), event_observed = logical(0))
  dat_E <- dat_T

  current <- cfg$start_dose
  now <- 0
  n_total <- 0
  path <- list()

  while (n_total < cfg$n_max) {
    # next patient's inter-arrival time (Poisson accrual), in days
    gap <- rexp(1, rate = cfg$accrual_rate / 7)
    now <- now + gap
    j <- current[1]; k <- current[2]

    # simulate this patient's true toxicity/efficacy event times (exponential for simplicity)
    has_tox <- rbinom(1, 1, true_pT[j, k])
    has_eff <- rbinom(1, 1, true_pE[j, k])
    t_time  <- if (has_tox) runif(1, 1, cfg$assessment_window_T) else NA
    e_time  <- if (has_eff) runif(1, 1, cfg$assessment_window_E) else NA

    dat_T <- rbind(dat_T, data.frame(dose1 = j, dose2 = k, enroll_time = now,
                                      event_time = t_time, event_observed = has_tox == 1))
    dat_E <- rbind(dat_E, data.frame(dose1 = j, dose2 = k, enroll_time = now,
                                      event_time = e_time, event_observed = has_eff == 1))
    n_total <- n_total + 1

    # decision uses only information matured "as of now" (right-censor at 'now')
    post_now <- utpi_comb_posterior(dat_T, dat_E, cfg, now = now)
    step <- utpi_comb_next_dose(current, post_now, cfg)
    path[[length(path) + 1]] <- list(dose = current, time = now, tox = has_tox, eff = has_eff)

    if (step$stop_trial) { message("Trial stopped: no admissible dose pairs remain."); break }
    current <- step$next_dose
  }

  # final analysis: let every patient fully mature (evaluate far in the future)
  final_time <- now + max(cfg$assessment_window_T, cfg$assessment_window_E) + 1
  final_post <- utpi_comb_posterior(dat_T, dat_E, cfg, now = final_time)
  A_final <- utpi_comb_admissible(final_post, cfg)

  if (any(A_final)) {
    idx <- which(A_final, arr.ind = TRUE)
    sel <- select_all_tied(idx, final_post$U_hat[A_final])   # ALL admissible cells tied on max U -- no tie-break
  } else {
    sel <- matrix(c(NA_real_, NA_real_), nrow = 1, dimnames = list(NULL, c("j", "k")))
  }

  list(dat_T = dat_T, dat_E = dat_E, posterior = final_post, admissible = A_final,
       recommended_OBDC = sel, n_tied_OBDC = nrow(sel), path = path, n_total = n_total)
}

## -----------------------------------------------------------------------------
## 5. Example
## -----------------------------------------------------------------------------

if (identical(environment(), globalenv())) {
  cfg <- utpi_comb_config(J = 4, K = 4, phiT = 0.30, phiE = 0.20, wE = 1, wT = 1,
                           accrual_rate = 2, n_max = 36)

  true_pT <- matrix(c(0.05, 0.10, 0.18, 0.30,
                       0.10, 0.18, 0.28, 0.40,
                       0.18, 0.28, 0.38, 0.50,
                       0.30, 0.40, 0.50, 0.62), nrow = 4, byrow = TRUE)
  true_pE <- matrix(c(0.10, 0.20, 0.30, 0.32,
                       0.20, 0.35, 0.45, 0.46,
                       0.30, 0.45, 0.55, 0.55,
                       0.32, 0.46, 0.55, 0.55), nrow = 4, byrow = TRUE)

  res <- simulate_utpi_comb(cfg, true_pE, true_pT, seed = 2026)
  if (anyNA(res$recommended_OBDC)) {
    cat("Recommended OBDC: none -- no admissible dose pair.\n")
  } else if (nrow(res$recommended_OBDC) == 1) {
    cat("Recommended OBDC (Agent A level, Agent B level):", res$recommended_OBDC[1, ], "\n")
  } else {
    cat("Recommended OBDC:", nrow(res$recommended_OBDC), "cells tied on utility (all reported, no tie-break):\n")
    print(res$recommended_OBDC)
  }
  cat("Total patients accrued:", res$n_total, "\n")
  cat("Posterior mean utility by dose pair:\n")
  print(round(res$posterior$U_hat, 3))
}
