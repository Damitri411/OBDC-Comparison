## =============================================================================
## Comb-BOIN12: model-assisted Bayesian design for two-agent OBDC selection
## Implements the design as formalised in Sec. 5.2 of Mukherjee, Takeda & Wason
## "Optimal Dose Combination Selection in Oncology Trials Using Bayesian Designs"
## (based on Lu, Zhang, Yuan & Lin, 2025, Stat Biopharm Res)
##
## NOTE: this is an original re-implementation built directly from the closed-form
## equations given in the manuscript. It is NOT the authors' validated software
## (no public code release for Comb-BOIN12 could be located). Treat as a design-
## exploration / teaching tool; validate independently before any operational use.
## =============================================================================

## -----------------------------------------------------------------------------
## 1. Design configuration
## -----------------------------------------------------------------------------

comb_boin12_config <- function(
  J = 4, K = 4,                 # dose grid: J levels of Agent A, K levels of Agent B
  phiT   = 0.30,                # target DLT rate
  phi1   = 0.60 * phiT,         # lower reference toxicity (highly likely underdosing)
  phi2   = 1.40 * phiT,         # upper reference toxicity (highly likely overdosing)
  phiE   = 0.20,                # efficacy floor
  wE     = 1, wT = 1,           # clinical utility weights, U = wE*pE - wT*pT
  aE = 1, bE = 1,               # Beta prior for efficacy, per cell
  aT = 1, bT = 1,               # Beta prior for toxicity, per cell
  alpha_u = 1, beta_u = 1,      # Beta prior for the scaled quasi-binomial utility
  CT = 0.90,                    # posterior cap: Pr(pT > phiT | data) <= CT -> admissible
  CE = 0.90,                    # posterior cap: Pr(pE < phiE | data) <= CE -> admissible
  cohort_size = 3,
  max_n_per_cell = 15,
  n_max = 60,                   # total sample size cap
  start_dose = c(1, 1)
) {
  lam <- boin_boundaries(phi1, phiT, phi2)
  list(J = J, K = K, phiT = phiT, phi1 = phi1, phi2 = phi2, phiE = phiE,
       wE = wE, wT = wT, aE = aE, bE = bE, aT = aT, bT = bT,
       alpha_u = alpha_u, beta_u = beta_u, CT = CT, CE = CE,
       lambda1 = lam$lambda1, lambda2 = lam$lambda2,
       cohort_size = cohort_size, max_n_per_cell = max_n_per_cell,
       n_max = n_max, start_dose = start_dose)
}

## Standard BOIN escalation/de-escalation boundaries (Yuan & Liu, 2015),
## used here for the toxicity-movement layer described in Sec. 5.2.
boin_boundaries <- function(phi1, phiT, phi2) {
  lambda1 <- log((1 - phi1) / (1 - phiT)) / log((phiT * (1 - phi1)) / (phi1 * (1 - phiT)))
  lambda2 <- log((1 - phiT) / (1 - phi2)) / log((phi2 * (1 - phiT)) / (phiT * (1 - phi2)))
  list(lambda1 = lambda1, lambda2 = lambda2)
}

## -----------------------------------------------------------------------------
## 2. Trial state: per-cell counts
## -----------------------------------------------------------------------------

new_trial_state <- function(cfg) {
  list(
    n   = matrix(0, cfg$J, cfg$K),
    nE  = matrix(0, cfg$J, cfg$K),   # efficacy responses
    nT  = matrix(0, cfg$J, cfg$K),   # DLTs
    tried = matrix(FALSE, cfg$J, cfg$K)
  )
}

## -----------------------------------------------------------------------------
## 3. Posterior summaries (Sec. 5.2 equations)
## -----------------------------------------------------------------------------

## Posterior mean/tail probabilities for pT, pE at every cell (independent
## Beta-Binomial conjugate updates), plus the quasi-binomial utility posterior:
##   x_jk = (nE_jk + wT*(n_jk - nT_jk)) / (1 + wT)
##   u*_jk = (U_jk + wT) / (1 + wT)  ~  Beta(alpha + x_jk, beta + n_jk - x_jk)
posterior_summary <- function(state, cfg) {
  J <- cfg$J; K <- cfg$K
  pT_hat <- pE_hat <- U_hat <- matrix(NA_real_, J, K)
  Pr_over_T <- Pr_under_E <- matrix(NA_real_, J, K)

  for (j in 1:J) for (k in 1:K) {
    n  <- state$n[j, k]; nE <- state$nE[j, k]; nT <- state$nT[j, k]

    pT_hat[j, k] <- (cfg$aT + nT) / (cfg$aT + cfg$bT + n)
    pE_hat[j, k] <- (cfg$aE + nE) / (cfg$aE + cfg$bE + n)

    Pr_over_T[j, k]  <- 1 - pbeta(cfg$phiT, cfg$aT + nT, cfg$bT + n - nT)
    Pr_under_E[j, k] <- pbeta(cfg$phiE, cfg$aE + nE, cfg$bE + n - nE)

    xjk <- (nE + cfg$wT * (n - nT)) / (1 + cfg$wT)
    a_post <- cfg$alpha_u + xjk
    b_post <- cfg$beta_u + n - xjk
    u_star_mean <- a_post / (a_post + b_post)          # posterior mean of scaled utility
    U_hat[j, k] <- u_star_mean * (1 + cfg$wT) - cfg$wT  # back-transform to U scale
  }

  list(pT_hat = pT_hat, pE_hat = pE_hat, U_hat = U_hat,
       Pr_over_T = Pr_over_T, Pr_under_E = Pr_under_E)
}

## Admissible set: A = { (j,k) : Pr(pT>phiT|D) <= CT  and  Pr(pE<phiE|D) <= CE }
admissible_set <- function(post, cfg, tried) {
  A <- (post$Pr_over_T <= cfg$CT) & (post$Pr_under_E <= cfg$CE)
  A[!tried] <- FALSE   # can only rank cells that have been tried in this simple version
  A
}

## -----------------------------------------------------------------------------
## 4. Local BOIN toxicity movement + utility ranking among neighbours
## -----------------------------------------------------------------------------

neighbours <- function(j, k, J, K) {
  cand <- rbind(c(j + 1, k), c(j, k + 1),   # escalate (one axis at a time: no-skip rule)
                c(j - 1, k), c(j, k - 1),   # de-escalate
                c(j, k))                    # stay
  cand <- cand[cand[, 1] >= 1 & cand[, 1] <= J & cand[, 2] >= 1 & cand[, 2] <= K, , drop = FALSE]
  cand
}

## Toxicity-driven local move at the current cell (j,k), following the BOIN rule
##   escalate if pT_hat < lambda1, de-escalate if pT_hat > lambda2, else stay
## then, among admissible neighbours reachable from that move, pick the one
## maximising E[U | data] (Sec. 5.2: "allocation favours the admissible neighbour
## maximising E[U_jk | D]").
next_dose_comb_boin12 <- function(current, state, cfg) {
  j <- current[1]; k <- current[2]
  post <- posterior_summary(state, cfg)
  A <- admissible_set(post, cfg, state$tried)

  pT_here <- post$pT_hat[j, k]
  if (state$n[j, k] == 0) pT_here <- cfg$phiT  # no data yet: neutral

  move <- if (pT_here < cfg$lambda1) "escalate" else if (pT_here > cfg$lambda2) "de-escalate" else "stay"

  cand <- neighbours(j, k, cfg$J, cfg$K)
  # restrict candidates to those consistent with the indicated move direction
  keep <- switch(move,
    "escalate"   = (cand[, 1] >= j & cand[, 2] >= k) & !(cand[, 1] == j & cand[, 2] == k),
    "de-escalate"= (cand[, 1] <= j & cand[, 2] <= k) & !(cand[, 1] == j & cand[, 2] == k),
    "stay"       = (cand[, 1] == j & cand[, 2] == k)
  )
  cand <- cand[keep, , drop = FALSE]
  if (nrow(cand) == 0) cand <- matrix(c(j, k), nrow = 1)

  # rank by expected utility among ADMISSIBLE cells; fall back to untried/neutral cell if none admissible
  util <- apply(cand, 1, function(z) post$U_hat[z[1], z[2]])
  is_adm <- apply(cand, 1, function(z) A[z[1], z[2]] || state$n[z[1], z[2]] == 0)
  cand_ok <- cand[is_adm, , drop = FALSE]
  util_ok <- util[is_adm]

  if (nrow(cand_ok) == 0) { cand_ok <- matrix(c(j, k), nrow = 1); util_ok <- post$U_hat[j, k] }
  best <- cand_ok[which.max(util_ok), ]
  list(next_dose = as.numeric(best), move = move, post = post, admissible = A)
}

## -----------------------------------------------------------------------------
## 5. Trial simulation
## -----------------------------------------------------------------------------

## true_pE, true_pT: J x K matrices of TRUE efficacy/toxicity probabilities (scenario)
simulate_comb_boin12 <- function(cfg, true_pE, true_pT, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  state <- new_trial_state(cfg)
  current <- cfg$start_dose
  total_n <- 0
  path <- list()

  while (total_n < cfg$n_max) {
    j <- current[1]; k <- current[2]
    if (state$n[j, k] >= cfg$max_n_per_cell) break

    m <- min(cfg$cohort_size, cfg$n_max - total_n)
    nE_new <- rbinom(1, m, true_pE[j, k])
    nT_new <- rbinom(1, m, true_pT[j, k])

    state$n[j, k]  <- state$n[j, k] + m
    state$nE[j, k] <- state$nE[j, k] + nE_new
    state$nT[j, k] <- state$nT[j, k] + nT_new
    state$tried[j, k] <- TRUE
    total_n <- total_n + m

    step <- next_dose_comb_boin12(current, state, cfg)
    path[[length(path) + 1]] <- list(dose = current, n = m, nE = nE_new, nT = nT_new, move = step$move)

    # global safety stop: if starting dose is itself over-toxic with high confidence
    post_check <- posterior_summary(state, cfg)
    if (post_check$Pr_over_T[1, 1] > cfg$CT && state$n[1, 1] >= cfg$cohort_size) {
      message("Trial stopped for safety: lowest dose combination is over-toxic.")
      break
    }

    current <- step$next_dose
  }

  final_post <- posterior_summary(state, cfg)
  A_final <- admissible_set(final_post, cfg, state$tried)
  if (any(A_final)) {
    idx <- which(A_final, arr.ind = TRUE)
    Uvals <- final_post$U_hat[A_final]
    sel <- idx[which.max(Uvals), ]
    names(sel) <- c("j", "k")
  } else {
    sel <- c(j = NA, k = NA)
  }

  list(state = state, posterior = final_post, admissible = A_final,
       recommended_OBDC = sel, path = path, n_total = total_n)
}

## -----------------------------------------------------------------------------
## 6. Example
## -----------------------------------------------------------------------------

if (identical(environment(), globalenv())) {
  cfg <- comb_boin12_config(J = 4, K = 4, phiT = 0.30, phiE = 0.20,
                             wE = 1, wT = 1, cohort_size = 3, n_max = 48)

  # toy true scenario: toxicity increases with both agents; efficacy plateaus
  true_pT <- matrix(c(0.05, 0.10, 0.18, 0.30,
                       0.10, 0.18, 0.28, 0.40,
                       0.18, 0.28, 0.38, 0.50,
                       0.30, 0.40, 0.50, 0.62), nrow = 4, byrow = TRUE)
  true_pE <- matrix(c(0.10, 0.20, 0.30, 0.32,
                       0.20, 0.35, 0.45, 0.46,
                       0.30, 0.45, 0.55, 0.55,
                       0.32, 0.46, 0.55, 0.55), nrow = 4, byrow = TRUE)

  res <- simulate_comb_boin12(cfg, true_pE, true_pT, seed = 2026)
  cat("Recommended OBDC (Agent A level, Agent B level):", res$recommended_OBDC, "\n")
  cat("Total patients used:", res$n_total, "\n")
  print(round(res$posterior$U_hat, 3))
}
