## =============================================================================
## COMIC (Combination Optimization in Multiple IndiCations): utility-integrated
## Bayesian design for two-agent OBDC selection across multiple indications.
## Implements the design as formalised in Sec. 5.3 of Mukherjee, Takeda & Wason
## "Optimal Dose Combination Selection in Oncology Trials Using Bayesian Designs"
## (based on Chen, Takeda & Yuan, 2025, Stat Med, e70107)
##
## NOTE: original re-implementation from the manuscript's equations; no public
## code release for COMIC could be located. Treat as a design-exploration tool.
## =============================================================================

## -----------------------------------------------------------------------------
## 1. Design configuration
## -----------------------------------------------------------------------------

## Six ordinal outcome categories (Sec. 5.3), best (1) to worst (6):
##  1: efficacious, no toxicity            4: non-efficacious, no toxicity
##  2: efficacious, tolerable low-grade tox 5: non-efficacious, tolerable tox
##  3: efficacious, DLT                     6: non-efficacious, DLT
comic_config <- function(
  J = 4, K = 4,
  psi = c(100, 80, 40, 30, 15, 0),  # clinician-elicited utility per category (Table 1)
  phi1_bar = 0.35, phi2_bar = 0.35, # toxicity safety thresholds for agents 1 & 2
  psi_bar_m = 40,                   # minimum acceptable utility for indication m
  Cs1 = 0.90, Cs2 = 0.90,           # posterior caps for toxicity exceedance
  Ce  = 0.90,                       # posterior cap for utility shortfall
  theta_b = 0.80, pd_cutoff = 0.80, # PD saturation threshold & posterior confidence cutoff
  dirichlet_prior = rep(1, 6),      # flat Dirichlet(1,...,1) prior over the 6 categories, per cell
  pd_prior = c(1, 1),               # Beta prior (a,b) for "PD saturated" indicator, per cell
  cohort_size = 3, n_max_stage1 = 30, n_max_stage2 = 24,
  start_dose = c(1, 1)
) {
  list(J = J, K = K, psi = psi, phi1_bar = phi1_bar, phi2_bar = phi2_bar,
       psi_bar_m = psi_bar_m, Cs1 = Cs1, Cs2 = Cs2, Ce = Ce,
       theta_b = theta_b, pd_cutoff = pd_cutoff,
       dirichlet_prior = dirichlet_prior, pd_prior = pd_prior,
       cohort_size = cohort_size, n_max_stage1 = n_max_stage1, n_max_stage2 = n_max_stage2,
       start_dose = start_dose)
}

## -----------------------------------------------------------------------------
## 2. Trial state (per indication): category counts, per-agent DLT counts, PD counts
## -----------------------------------------------------------------------------

new_comic_state <- function(cfg) {
  list(
    n      = matrix(0, cfg$J, cfg$K),
    counts = array(0, dim = c(cfg$J, cfg$K, 6)),  # counts in each of the 6 categories
    n_dlt1 = matrix(0, cfg$J, cfg$K),             # DLT attributable to agent-1 axis (for phi1_bar check)
    n_dlt2 = matrix(0, cfg$J, cfg$K),             # DLT attributable to agent-2 axis (for phi2_bar check)
    n_pd_sat = matrix(0, cfg$J, cfg$K)            # number of patients with PD saturation observed
  )
}

## -----------------------------------------------------------------------------
## 3. Posterior updating (Dirichlet-multinomial conjugacy, Sec. 5.3)
##    mu*_jk = sum_l psi_l * pi_jk(l),   x_jk ~ Binomial(n_jk, mu_jk),
##    mu_jk | x_jk ~ Beta(1+x_jk, 1+n_jk-x_jk),  mu_jk = mu*_jk/100
## The manuscript gives both a full Dirichlet-multinomial route (category
## probabilities pi_jk(l)) and a quasi-binomial shortcut for the *mean* utility;
## we implement the Dirichlet-multinomial version (richer: also gives per-category
## posteriors, needed for the toxicity/PD admissibility checks) and use it to
## reconstruct mu*_jk = sum psi_l * E[pi_jk(l) | data].
## -----------------------------------------------------------------------------

comic_posterior_summary <- function(state, cfg) {
  J <- cfg$J; K <- cfg$K
  mu_hat <- matrix(NA_real_, J, K)          # posterior mean utility (0-100 scale)
  PrTox1_high <- PrTox2_high <- matrix(NA_real_, J, K)
  PrUtil_low  <- matrix(NA_real_, J, K)
  PD_sat_prob <- matrix(NA_real_, J, K)

  for (j in 1:J) for (k in 1:K) {
    a_post <- cfg$dirichlet_prior + state$counts[j, k, ]
    pi_hat <- a_post / sum(a_post)                    # posterior mean category probabilities
    mu_hat[j, k] <- sum(cfg$psi * pi_hat)

    n_jk <- state$n[j, k]
    # posterior for "toxicity on agent-1 axis rate" and "agent-2 axis rate" via Beta-Binomial
    PrTox1_high[j, k] <- 1 - pbeta(cfg$phi1_bar, 1 + state$n_dlt1[j, k], 1 + n_jk - state$n_dlt1[j, k])
    PrTox2_high[j, k] <- 1 - pbeta(cfg$phi2_bar, 1 + state$n_dlt2[j, k], 1 + n_jk - state$n_dlt2[j, k])

    # posterior for Pr(utility < psi_bar_m | data): approximate via Monte Carlo from
    # the Dirichlet posterior (closed form is not available for a linear combination
    # of Dirichlet components in general)
    if (n_jk > 0 || any(a_post != cfg$dirichlet_prior)) {
      draws <- MCMCpack_rdirichlet(2000, a_post)
      util_draws <- draws %*% cfg$psi
      PrUtil_low[j, k] <- mean(util_draws < cfg$psi_bar_m)
    } else {
      PrUtil_low[j, k] <- 0.5   # no information yet: neutral
    }

    PD_sat_prob[j, k] <- 1 - pbeta(cfg$theta_b, cfg$pd_prior[1] + state$n_pd_sat[j, k],
                                    cfg$pd_prior[2] + n_jk - state$n_pd_sat[j, k])
  }

  list(mu_hat = mu_hat, PrTox1_high = PrTox1_high, PrTox2_high = PrTox2_high,
       PrUtil_low = PrUtil_low, PD_sat_prob = PD_sat_prob)
}

## Minimal, dependency-free Dirichlet sampler (equivalent to MCMCpack::rdirichlet)
MCMCpack_rdirichlet <- function(n, alpha) {
  k <- length(alpha)
  x <- matrix(rgamma(n * k, shape = alpha, rate = 1), nrow = n, byrow = TRUE)
  x / rowSums(x)
}

## Admissible set (Sec. 5.3):
##   Pr(tox1 > phi1_bar) <= Cs1,  Pr(tox2 > phi2_bar) <= Cs2,  Pr(utility < psi_bar_m) <= Ce
comic_admissible_set <- function(post, cfg) {
  (post$PrTox1_high <= cfg$Cs1) & (post$PrTox2_high <= cfg$Cs2) & (post$PrUtil_low <= cfg$Ce)
}

## -----------------------------------------------------------------------------
## 4. PD-guided axis-wise escalation
##    Once Pr(PD saturated | data) crosses pd_cutoff on one agent's axis, that
##    agent holds while the other escalates (Sec. 5.3).
## -----------------------------------------------------------------------------

comic_neighbours <- function(j, k, J, K, agent1_saturated, agent2_saturated) {
  cand <- list()
  if (!agent1_saturated && j < J) cand[[length(cand) + 1]] <- c(j + 1, k)
  if (!agent2_saturated && k < K) cand[[length(cand) + 1]] <- c(j, k + 1)
  if (j > 1) cand[[length(cand) + 1]] <- c(j - 1, k)
  if (k > 1) cand[[length(cand) + 1]] <- c(j, k - 1)
  cand[[length(cand) + 1]] <- c(j, k)
  do.call(rbind, cand)
}

comic_next_dose <- function(current, state, cfg) {
  j <- current[1]; k <- current[2]
  post <- comic_posterior_summary(state, cfg)
  A <- comic_admissible_set(post, cfg)

  agent1_saturated <- post$PD_sat_prob[j, k] >= cfg$pd_cutoff  # agent-1 axis "full", hold agent 1
  agent2_saturated <- FALSE  # symmetric extension left to the user's PD assay design; see README

  cand <- comic_neighbours(j, k, cfg$J, cfg$K, agent1_saturated, agent2_saturated)
  is_adm <- apply(cand, 1, function(z) A[z[1], z[2]] || state$n[z[1], z[2]] == 0)
  cand_ok <- cand[is_adm, , drop = FALSE]
  if (nrow(cand_ok) == 0) { cand_ok <- matrix(c(j, k), nrow = 1) }

  util <- apply(cand_ok, 1, function(z) post$mu_hat[z[1], z[2]])
  best <- cand_ok[which.max(util), ]
  list(next_dose = as.numeric(best), post = post, admissible = A)
}

## -----------------------------------------------------------------------------
## 5. Two-stage, multi-indication simulation
##    Stage I: escalate for ONE indication starting at low doses.
##    Stage II: initialise subsequent indications at the Stage-I provisional OBDC,
##    borrowing the Stage-I Dirichlet posterior as an informative prior (a simple,
##    transparent form of the "information borrowing" described in Sec. 5.3).
## -----------------------------------------------------------------------------

comic_run_indication <- function(cfg, true_cat_probs, true_dlt1, true_dlt2, true_pd_sat,
                                  n_max, start_dose, prior_override = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  state <- new_comic_state(cfg)
  if (!is.null(prior_override)) cfg$dirichlet_prior <- prior_override

  current <- start_dose; total_n <- 0; path <- list()
  while (total_n < n_max) {
    j <- current[1]; k <- current[2]
    m <- min(cfg$cohort_size, n_max - total_n)

    cats <- sample(1:6, m, replace = TRUE, prob = true_cat_probs[j, k, ])
    for (l in 1:6) state$counts[j, k, l] <- state$counts[j, k, l] + sum(cats == l)
    dlt1 <- rbinom(1, m, true_dlt1[j, k]); dlt2 <- rbinom(1, m, true_dlt2[j, k])
    pdsat <- rbinom(1, m, true_pd_sat[j, k])

    state$n[j, k] <- state$n[j, k] + m
    state$n_dlt1[j, k] <- state$n_dlt1[j, k] + dlt1
    state$n_dlt2[j, k] <- state$n_dlt2[j, k] + dlt2
    state$n_pd_sat[j, k] <- state$n_pd_sat[j, k] + pdsat
    total_n <- total_n + m

    step <- comic_next_dose(current, state, cfg)
    path[[length(path) + 1]] <- list(dose = current, n = m, cats = table(factor(cats, 1:6)))
    current <- step$next_dose
  }

  post <- comic_posterior_summary(state, cfg)
  A <- comic_admissible_set(post, cfg)
  if (any(A)) {
    idx <- which(A, arr.ind = TRUE)
    sel <- idx[which.max(post$mu_hat[A]), ]; names(sel) <- c("j", "k")
  } else sel <- c(j = NA, k = NA)

  list(state = state, posterior = post, admissible = A, recommended_OBDC = sel, path = path)
}

## -----------------------------------------------------------------------------
## 6. Example: two indications, Stage I then Stage II with borrowing
## -----------------------------------------------------------------------------

if (identical(environment(), globalenv())) {
  cfg <- comic_config(J = 4, K = 4, cohort_size = 3, n_max_stage1 = 24, n_max_stage2 = 18)

  make_true_cats <- function(pE, pT_dlt, pT_lowgrade) {
    # collapse (pE,pT) scenario into the 6-category probability array
    J <- nrow(pE); K <- ncol(pE)
    arr <- array(0, dim = c(J, K, 6))
    for (j in 1:J) for (k in 1:K) {
      pE_jk <- pE[j, k]; pDLT <- pT_dlt[j, k]; pLow <- pT_lowgrade[j, k]
      pNoTox <- max(0, 1 - pDLT - pLow)
      arr[j, k, ] <- c(pE_jk * pNoTox, pE_jk * pLow, pE_jk * pDLT,
                        (1 - pE_jk) * pNoTox, (1 - pE_jk) * pLow, (1 - pE_jk) * pDLT)
      arr[j, k, ] <- arr[j, k, ] / sum(arr[j, k, ])
    }
    arr
  }

  pE  <- matrix(c(0.10,0.25,0.40,0.45, 0.20,0.40,0.55,0.58, 0.30,0.50,0.60,0.60, 0.32,0.52,0.60,0.60), 4, 4, byrow=TRUE)
  pDLT<- matrix(c(0.03,0.06,0.10,0.18, 0.06,0.10,0.18,0.28, 0.10,0.18,0.28,0.40, 0.18,0.28,0.40,0.55), 4, 4, byrow=TRUE)
  pLow<- matrix(0.10, 4, 4)
  true_cats <- make_true_cats(pE, pDLT, pLow)
  true_dlt1 <- pDLT; true_dlt2 <- pDLT           # simplified: shared DLT surface across both agent axes
  true_pd_sat <- matrix(c(0.05,0.15,0.40,0.75, 0.10,0.30,0.65,0.90, 0.20,0.50,0.85,0.95, 0.30,0.60,0.90,0.97), 4, 4, byrow=TRUE)

  cat("=== Stage I: indication 1 ===\n")
  res1 <- comic_run_indication(cfg, true_cats, true_dlt1, true_dlt2, true_pd_sat,
                                n_max = cfg$n_max_stage1, start_dose = cfg$start_dose, seed = 11)
  cat("Stage-I provisional OBDC:", res1$recommended_OBDC, "\n")
  print(round(res1$posterior$mu_hat, 1))

  cat("\n=== Stage II: indication 2 (borrows Stage-I posterior, starts at provisional OBDC) ===\n")
  prior_borrow <- cfg$dirichlet_prior + apply(res1$state$counts, 3, sum) / cfg$J / cfg$K  # down-weighted borrowing
  start2 <- if (!any(is.na(res1$recommended_OBDC))) res1$recommended_OBDC else cfg$start_dose
  res2 <- comic_run_indication(cfg, true_cats, true_dlt1, true_dlt2, true_pd_sat,
                                n_max = cfg$n_max_stage2, start_dose = start2,
                                prior_override = prior_borrow, seed = 22)
  cat("Stage-II (indication 2) OBDC:", res2$recommended_OBDC, "\n")
  print(round(res2$posterior$mu_hat, 1))
}
