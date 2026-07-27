## =============================================================================
## execution.R  --  RUN THIS FILE:   Rscript execution.R
## =============================================================================
## A single, self-contained script that runs and compares five OBDC designs
## (Comb-BOIN12, EffTox, EffTox-TITE, COMIC, uTPI-Comb) from one common
## parameter set, across one or more dose-combination scenarios, and prints a
## report for every scenario x design combination.
##
## EffTox-TITE is the same EffTox model with continuous (rolling) accrual and
## a time-to-event weighted likelihood (isTITE = TRUE in efftox_config()),
## instead of EffTox's usual fixed-cohort, fully-observed-outcome conduct.
##
## For every scenario, the report also prints the TRUE admissible set and
## TRUE OBDC(s) from ground truth -- if multiple cells tie exactly on utility
## among admissible cells, ALL of them are reported, never silently reduced
## to one "winner". The same no-tie-break rule applies inside each design's
## own simulated trial (see efftox.R etc.: a trial's final recommended_OBDC
## is now a full tied set, not one cell) -- the Operating Characteristics
## report below aggregates many such trials, so it additionally shows a
## "modal cell across replicates" as a display summary (which IS tie-broken,
## purely to have one cell to highlight) alongside how often an individual
## trial's own final selection was itself untied vs. tied.
##
## Built around Table 2 of Mukherjee, Takeda & Wason (manuscript) and the toy
## scenarios from Table 1 of the uTPI-Comb paper (Liang, Yang & Yuan, 2024,
## JUSTC, DOI: 10.52396/JUSTC-2024-0104). Full references are in README.md.
## An interactive version of the same comparisons is available online -- see
## the link in README.md.
##
## HOW TO USE THIS FILE
##   1. Edit the "SETTINGS YOU CAN CHANGE" block just below.
##   2. Save, then run from a terminal in this folder:   Rscript execution.R
##   3. Read the report that prints to the console (or redirect to a file:
##        Rscript execution.R > my_report.txt
##      ).
##
## It sources four sibling files that must sit in the same folder:
##   comb_boin12.R, efftox.R, comic.R, utpi_comb.R   (the four design engines --
##     efftox.R contains BOTH the usual EffTox and the EffTox-TITE variant)
##   scenarios.R                                      (the five fixed scenarios)
## =============================================================================


## #############################################################################
## SETTINGS YOU CAN CHANGE
## #############################################################################

## Which of the five fixed uTPI-Comb-paper scenarios to run (comment lines out
## to skip a scenario, or add "Custom" -- see CUSTOM_SCENARIO below):
SCENARIOS_TO_RUN <- c(
  "Scenario 1",
  "Scenario 3",
  "Scenario 5",
  "Scenario 9",
  "Scenario 10"
  # , "Custom"          # <- uncomment to also run the custom scenario below
)

## Table 2 parameters (Sec. 6.2) -- the common settings every design and every
## scenario uses:
DELTA_E <- 0.20     # Minimal efficacy: threshold for acceptable efficacy
DELTA_T <- 0.30     # Maximal toxicity: threshold for acceptable toxicity
PI_E    <- 0.80     # Posterior requirement (efficacy): Pr(E >= deltaE | data) >= piE
PI_T    <- 0.80     # Posterior requirement (toxicity): Pr(T <= deltaT | data) >= piT
W_E     <- 1        # Utility weight on efficacy
W_T     <- 1        # Utility weight on toxicity
A_E     <- 1; B_E <- 1   # Efficacy Beta prior (aE, bE)
A_T     <- 1; B_T <- 1   # Toxicity Beta prior (aT, bT)
N_PER_CELL <- 6     # OC per-cell sample size (patients per dose-grid cell)
N_SIM      <- 10    # OC replicates (more = slower but more precise)

## Optional: define your own scenario (only used if "Custom" is included in
## SCENARIOS_TO_RUN above). J,K must each be >= 2. Leave true_pT/true_pE as
## NULL to get an editable smooth default; or fill in your own J x K matrices.
CUSTOM_J <- 4
CUSTOM_K <- 4
CUSTOM_TRUE_PT <- NULL   # e.g. matrix(c(0.05,0.10,..., ...), nrow=CUSTOM_J, byrow=TRUE)
CUSTOM_TRUE_PE <- NULL

## #############################################################################
## Nothing below this line needs editing for routine use.
## #############################################################################

SRC_DIR <- getwd()

## Friendly pre-flight check: if you're new to R, the single most common
## reason this script fails is that R's working directory isn't set to the
## folder containing these files. This checks for that up front and tells
## you exactly what's missing and how to fix it, instead of a cryptic
## "cannot open file" error partway through.
required_siblings <- c("scenarios.R", "comb_boin12.R", "efftox.R", "comic.R", "utpi_comb.R")
missing_siblings <- required_siblings[!file.exists(file.path(SRC_DIR, required_siblings))]
if (length(missing_siblings) > 0) {
  stop(
    "\n\nCan't find these required file(s) in the current working directory:\n  ",
    paste(missing_siblings, collapse = "\n  "),
    "\n\nCurrent working directory is:\n  ", SRC_DIR,
    "\n\nFix: put execution.R and all its sibling files (scenarios.R, comb_boin12.R,",
    "\nefftox.R, comic.R, utpi_comb.R) in the SAME folder, then either:",
    "\n  - in RStudio: Session -> Set Working Directory -> Choose Directory... (pick that folder), or",
    "\n  - in R:       setwd(\"path/to/that/folder\")",
    "\nand run execution.R again.\n",
    call. = FALSE
  )
}

source(file.path(SRC_DIR, "scenarios.R"))        # -> UTPI_PAPER_SCENARIOS, make_custom_scenario()

env_boin  <- new.env()
env_et    <- new.env()
env_comic <- new.env()
env_utpi  <- new.env()
sys.source(file.path(SRC_DIR, "comb_boin12.R"), envir = env_boin)
sys.source(file.path(SRC_DIR, "efftox.R"),      envir = env_et)
sys.source(file.path(SRC_DIR, "comic.R"),       envir = env_comic)
sys.source(file.path(SRC_DIR, "utpi_comb.R"),   envir = env_utpi)

CFG <- list(deltaE = DELTA_E, deltaT = DELTA_T, piE = PI_E, piT = PI_T,
            wE = W_E, wT = W_T, aE = A_E, bE = B_E, aT = A_T, bT = B_T,
            n = N_PER_CELL, Nsim = N_SIM)

## Assemble the requested scenarios into one lookup list of {J,K,true_pT,true_pE,source}
ALL_SCENARIOS <- UTPI_PAPER_SCENARIOS
ALL_SCENARIOS[["Custom"]] <- make_custom_scenario(CUSTOM_J, CUSTOM_K, CUSTOM_TRUE_PT, CUSTOM_TRUE_PE)

## =============================================================================
## Design note on piE/piT (worth reading once):
## Table 2's (piE,piT) are a FINAL-ANALYSIS admissible-set criterion (Sec. 6.3):
## Pr(E>=deltaE|data)>=piE and Pr(T<=deltaT|data)>=piT, evaluated once at the
## end of a trial. Each design's own SEQUENTIAL conduct (Sec. 5) has its own,
## separately-tuned interim safety/futility caps, which we leave at each
## design's published/script defaults so a trial can actually accrue past its
## first cohort. piE/piT are applied below purely to compute the reported
## "membership in A" and "mean |A|" at the final look.
## =============================================================================

build_cfg_comb_boin12 <- function(mc) {
  env_boin$comb_boin12_config(
    J = mc$J, K = mc$K, phiT = mc$deltaT, phiE = mc$deltaE,
    wE = mc$wE, wT = mc$wT, aE = mc$aE, bE = mc$bE, aT = mc$aT, bT = mc$bT,
    cohort_size = 3, n_max = mc$n * mc$J * mc$K, start_dose = c(1, 1)
  )
}
build_cfg_efftox <- function(mc, fast = TRUE) {
  env_et$efftox_config(
    dose1 = seq_len(mc$J), dose2 = seq_len(mc$K), thetaE = mc$deltaE, thetaT = mc$deltaT,
    alphaT = 0.70, alphaE = 0.70, n_iter = if (fast) 1200 else 4000, n_burnin = if (fast) 400 else 1500
  )
}
build_cfg_efftox_tite <- function(mc, fast = TRUE) {
  ## Same EffTox model as build_cfg_efftox(), but with isTITE = TRUE: continuous
  ## (Poisson) accrual and a time-to-event weighted likelihood instead of fixed
  ## cohorts of fully-observed outcomes. Tmax/accrual_rate are illustrative
  ## defaults (a 42-unit assessment window, 2 patients per unit time to match
  ## uTPI-Comb's accrual_rate elsewhere in this file) -- edit as needed for a
  ## real trial's actual assessment window and expected accrual speed.
  env_et$efftox_config(
    dose1 = seq_len(mc$J), dose2 = seq_len(mc$K), thetaE = mc$deltaE, thetaT = mc$deltaT,
    alphaT = 0.70, alphaE = 0.70, n_iter = if (fast) 1200 else 4000, n_burnin = if (fast) 400 else 1500,
    isTITE = TRUE, Tmax = 42, accrual_rate = 2
  )
}
build_cfg_comic <- function(mc) {
  psi_bar_m <- 100 * max(0, mc$wE * mc$deltaE - mc$wT * mc$deltaT)
  env_comic$comic_config(
    J = mc$J, K = mc$K, phi1_bar = mc$deltaT, phi2_bar = mc$deltaT,
    psi_bar_m = psi_bar_m, pd_prior = c(mc$aT, mc$bT), cohort_size = 3,
    n_max_stage1 = mc$n * mc$J * mc$K, n_max_stage2 = round(mc$n * mc$J * mc$K * 0.7), start_dose = c(1, 1)
  )
}
build_cfg_utpi_comb <- function(mc) {
  env_utpi$utpi_comb_config(
    J = mc$J, K = mc$K, phiT = mc$deltaT, phiE = mc$deltaE,
    wE = mc$wE, wT = mc$wT, aE = mc$aE, bE = mc$bE, aT = mc$aT, bT = mc$bT,
    mu_b_source = c(phiE = mc$deltaE, phiT = mc$deltaT), accrual_rate = 2, n_max = mc$n * mc$J * mc$K
  )
}

table2_admissible_comb_boin12 <- function(post, mc) (1 - post$Pr_under_E >= mc$piE) & (1 - post$Pr_over_T >= mc$piT)
table2_admissible_efftox      <- function(post, mc) (post$PrE_ok >= mc$piE) & (post$PrT_ok >= mc$piT)
table2_admissible_comic       <- function(post, mc) (1 - post$PrUtil_low >= mc$piE) & (1 - post$PrTox1_high >= mc$piT)
table2_admissible_utpi_comb   <- function(post, mc) (post$PrE_ok >= mc$piE) & (1 - post$PrT_over >= mc$piT)

## Multiple cells CAN legitimately tie. Two different situations use two
## different policies:
##  - The "modal cell across Nsim replicates" below is a SELECTION-FREQUENCY
##    tie (how often a cell was part of the recommended set across many
##    simulated trials) -- a summary display statistic, not a trial decision,
##    but it still needs to resolve to one highlighted cell for the printed
##    report, so it uses select_best_with_tiebreak() (lowest total dose j+k)
##    and reports how many cells were tied.
##  - The TRUE OBDC (ground truth) and each design's own FINAL recommended
##    OBDC (computed inside the design engines) are REPORTED CONCLUSIONS, not
##    operational necessities -- so neither is tie-broken. Every admissible
##    cell achieving the max utility is reported, via select_all_tied()
##    (defined in every engine file too) or its use below.
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

## =============================================================================
## TRUE admissible set and TRUE OBDC(s) -- computed directly from a scenario's
## ground-truth surfaces (no posterior/simulation involved). Reports EVERY
## admissible cell tied for the highest true utility U = wE*true_pE -
## wT*true_pT; no tie-break, since this is a factual property of the
## scenario itself, not an operational trial decision.
## =============================================================================

true_admissible_and_obdc <- function(mc, true_pE, true_pT) {
  A_true <- (true_pE >= mc$deltaE) & (true_pT <= mc$deltaT)
  U_true <- mc$wE * true_pE - mc$wT * true_pT
  if (!any(A_true)) {
    return(list(A_true = A_true, U_true = U_true,
                obdc = matrix(NA_real_, 0, 2, dimnames = list(NULL, c("j", "k")))))
  }
  idx <- which(A_true, arr.ind = TRUE)
  u_at_idx <- U_true[A_true]
  obdc <- select_all_tied(idx, u_at_idx)
  list(A_true = A_true, U_true = U_true, obdc = obdc)
}

## Prints the true admissible set and ALL true-OBDC candidates for one
## scenario (ground truth -- independent of any design).
report_true_obdc <- function(scenario_label, scenario, mc) {
  ta <- true_admissible_and_obdc(mc, scenario$true_pE, scenario$true_pT)
  cat("\n--- Ground truth for", scenario_label, "---\n")
  A_cells <- which(ta$A_true, arr.ind = TRUE)
  if (nrow(A_cells) == 0) {
    cat("True admissible set A = {}  (empty) -- no cell satisfies true_pE >= deltaE and true_pT <= deltaT.\n\n")
    return(invisible(NULL))
  }
  cells_str <- paste0("(", A_cells[, "row"], ",", A_cells[, "col"], ")", collapse = ", ")
  cat("True admissible set A = {", cells_str, "}\n")
  if (nrow(ta$obdc) == 1) {
    j <- ta$obdc[1, "j"]; k <- ta$obdc[1, "k"]
    cat(sprintf("True OBDC = (j=%d, k=%d)   true_pE=%.3f  true_pT=%.3f  U=%.3f\n",
                j, k, scenario$true_pE[j, k], scenario$true_pT[j, k], ta$U_true[j, k]))
  } else {
    cat("True OBDC:", nrow(ta$obdc), "cells tied on true utility U (ALL reported, no tie-break):\n")
    for (r in seq_len(nrow(ta$obdc))) {
      j <- ta$obdc[r, "j"]; k <- ta$obdc[r, "k"]
      cat(sprintf("  candidate %d/%d: (j=%d, k=%d)   true_pE=%.3f  true_pT=%.3f  U=%.3f\n",
                  r, nrow(ta$obdc), j, k, scenario$true_pE[j, k], scenario$true_pT[j, k], ta$U_true[j, k]))
    }
  }
  cat("\n")
}

## `reps` is a list where each element's `$sel` is now a MATRIX (>=1 rows) of
## every cell tied for that replicate's own final recommendation -- since the
## design engines no longer tie-break their final OBDC. Every tied cell
## contributes to the heatmap and to the Ehat/That/Uhat/etc. summary stats
## below (so a replicate with 2 tied cells contributes 2 "cell instances",
## not 1) -- this is the most transparent way to aggregate when the
## underlying unit (a single trial's answer) can itself be a set.
summarise_oc <- function(reps, J, K, true_pE, true_pT, deltaE, deltaT) {
  Nsim <- length(reps); heat <- matrix(0, J, K); no_selection <- 0
  Ehat_sel <- That_sel <- Uhat_sel <- Etrue_sel <- Ttrue_sel <- A_sizes <- numeric(0)
  member_flags <- logical(0); n_unsafe <- 0; n_ineffective <- 0
  replicate_tie_flags <- logical(0)
  n_cell_instances <- 0
  
  for (r in reps) {
    A_sizes <- c(A_sizes, r$A_size)
    sel <- r$sel
    if (is.null(sel) || nrow(sel) == 0 || anyNA(sel)) { no_selection <- no_selection + 1; next }
    n_tied_this_rep <- nrow(sel)
    replicate_tie_flags <- c(replicate_tie_flags, n_tied_this_rep > 1)
    n_cell_instances <- n_cell_instances + n_tied_this_rep
    for (rr in seq_len(n_tied_this_rep)) heat[sel[rr, 1], sel[rr, 2]] <- heat[sel[rr, 1], sel[rr, 2]] + 1
    Ehat_sel <- c(Ehat_sel, r$Ehat); That_sel <- c(That_sel, r$That); Uhat_sel <- c(Uhat_sel, r$Uhat)
    Etrue_sel <- c(Etrue_sel, true_pE[sel]); Ttrue_sel <- c(Ttrue_sel, true_pT[sel])
    member_flags <- c(member_flags, r$table2_member)
    n_unsafe <- n_unsafe + sum(true_pT[sel] > deltaT)
    n_ineffective <- n_ineffective + sum(true_pE[sel] < deltaE)
  }
  heat_prop <- heat / Nsim
  n_selected <- Nsim - no_selection
  n_tied_modal <- 0
  if (n_selected > 0 && max(heat) > 0) {
    idx <- which(heat == max(heat), arr.ind = TRUE)
    pick <- select_best_with_tiebreak(idx, rep(1, nrow(idx)))   # all tied cells already share max(heat)
    modal_cell <- pick$sel
    n_tied_modal <- nrow(idx)
  } else {
    modal_cell <- c(NA, NA)
  }
  list(heatmap = heat_prop, prop_no_selection = no_selection / Nsim, modal_cell = modal_cell,
       n_tied_modal = n_tied_modal,
       prop_replicates_with_tie = if (n_selected > 0) mean(replicate_tie_flags) else NA,
       Ehat_at_modal = if (n_cell_instances > 0) mean(Ehat_sel) else NA,
       That_at_modal = if (n_cell_instances > 0) mean(That_sel) else NA,
       Uhat_at_modal = if (n_cell_instances > 0) mean(Uhat_sel) else NA,
       prop_modal_table2_member = if (n_cell_instances > 0) mean(member_flags) else NA,
       Pr_select_unsafe = if (n_cell_instances > 0) n_unsafe / n_cell_instances else NA,
       Pr_select_ineffective = if (n_cell_instances > 0) n_ineffective / n_cell_instances else NA,
       mean_true_E_at_selection = if (n_cell_instances > 0) mean(Etrue_sel) else NA,
       mean_true_T_at_selection = if (n_cell_instances > 0) mean(Ttrue_sel) else NA,
       mean_A_size = if (length(A_sizes) > 0) mean(A_sizes) else NA)
}

run_oc_comb_boin12 <- function(mc, true_pE, true_pT, seeds) {
  cfg <- build_cfg_comb_boin12(mc)
  reps <- lapply(seeds, function(s) {
    res <- env_boin$simulate_comb_boin12(cfg, true_pE, true_pT, seed = s)
    A_t2 <- table2_admissible_comb_boin12(res$posterior, mc)
    sel <- res$recommended_OBDC   # matrix: every cell tied for this trial's own final recommendation
    if (anyNA(sel)) return(list(sel = sel, A_size = sum(A_t2)))
    list(sel = sel, Ehat = res$posterior$pE_hat[sel], That = res$posterior$pT_hat[sel],
         Uhat = res$posterior$U_hat[sel], table2_member = A_t2[sel], A_size = sum(A_t2))
  })
  summarise_oc(reps, mc$J, mc$K, true_pE, true_pT, mc$deltaE, mc$deltaT)
}
run_oc_efftox <- function(mc, true_pE, true_pT, seeds, fast = TRUE) {
  cfg <- build_cfg_efftox(mc, fast = fast)
  reps <- lapply(seeds, function(s) {
    res <- env_et$simulate_efftox(cfg, true_pE, true_pT, cohort_size = 3, n_max = mc$n * mc$J * mc$K, seed = s)
    A_t2 <- table2_admissible_efftox(res$posterior, mc)
    sel <- res$recommended_OBDC
    if (anyNA(sel)) return(list(sel = sel, A_size = sum(A_t2)))
    list(sel = sel, Ehat = res$posterior$pE_hat[sel], That = res$posterior$pT_hat[sel],
         Uhat = res$posterior$Dbar[sel], table2_member = A_t2[sel], A_size = sum(A_t2))
  })
  summarise_oc(reps, mc$J, mc$K, true_pE, true_pT, mc$deltaE, mc$deltaT)
}
run_oc_efftox_tite <- function(mc, true_pE, true_pT, seeds, fast = TRUE) {
  ## Same as run_oc_efftox(), but with the isTITE = TRUE config: continuous
  ## accrual and TITE-weighted interim decisions. simulate_efftox() dispatches
  ## to the TITE path internally based on cfg$isTITE, so this is otherwise
  ## identical -- same posterior fields (PrE_ok/PrT_ok/Dbar/etc.), so it can
  ## reuse table2_admissible_efftox() unchanged. cohort_size is accepted but
  ## ignored on the TITE path (accrual is one patient at a time, continuously).
  cfg <- build_cfg_efftox_tite(mc, fast = fast)
  reps <- lapply(seeds, function(s) {
    res <- env_et$simulate_efftox(cfg, true_pE, true_pT, n_max = mc$n * mc$J * mc$K, seed = s)
    A_t2 <- table2_admissible_efftox(res$posterior, mc)
    sel <- res$recommended_OBDC
    if (anyNA(sel)) return(list(sel = sel, A_size = sum(A_t2)))
    list(sel = sel, Ehat = res$posterior$pE_hat[sel], That = res$posterior$pT_hat[sel],
         Uhat = res$posterior$Dbar[sel], table2_member = A_t2[sel], A_size = sum(A_t2))
  })
  summarise_oc(reps, mc$J, mc$K, true_pE, true_pT, mc$deltaE, mc$deltaT)
}
run_oc_comic <- function(mc, true_pE, true_pT, seeds) {
  cfg <- build_cfg_comic(mc)
  J <- mc$J; K <- mc$K
  arr <- array(0, dim = c(J, K, 6))
  for (j in 1:J) for (k in 1:K) {
    pE_jk <- true_pE[j, k]; pDLT <- true_pT[j, k]; pLow <- 0.10
    pNoTox <- max(0, 1 - pDLT - pLow)
    arr[j, k, ] <- c(pE_jk * pNoTox, pE_jk * pLow, pE_jk * pDLT,
                     (1 - pE_jk) * pNoTox, (1 - pE_jk) * pLow, (1 - pE_jk) * pDLT)
    arr[j, k, ] <- arr[j, k, ] / sum(arr[j, k, ])
  }
  true_pd_sat <- matrix(pmin(0.97, true_pT + 0.3), J, K)
  reps <- lapply(seeds, function(s) {
    res <- env_comic$comic_run_indication(cfg, arr, true_pT, true_pT, true_pd_sat,
                                          n_max = cfg$n_max_stage1, start_dose = cfg$start_dose, seed = s)
    A_t2 <- table2_admissible_comic(res$posterior, mc)
    sel <- res$recommended_OBDC
    if (anyNA(sel)) return(list(sel = sel, A_size = sum(A_t2)))
    mu <- res$posterior$mu_hat[sel] / 100
    list(sel = sel, Ehat = true_pE[sel], That = true_pT[sel],
         Uhat = mu, table2_member = A_t2[sel], A_size = sum(A_t2))
  })
  summarise_oc(reps, mc$J, mc$K, true_pE, true_pT, mc$deltaE, mc$deltaT)
}
run_oc_utpi_comb <- function(mc, true_pE, true_pT, seeds) {
  cfg <- build_cfg_utpi_comb(mc)
  reps <- lapply(seeds, function(s) {
    res <- env_utpi$simulate_utpi_comb(cfg, true_pE, true_pT, seed = s)
    A_t2 <- table2_admissible_utpi_comb(res$posterior, mc)
    sel <- res$recommended_OBDC
    if (anyNA(sel)) return(list(sel = sel, A_size = sum(A_t2)))
    list(sel = sel, Ehat = res$posterior$pE_hat[sel], That = res$posterior$pT_hat[sel],
         Uhat = res$posterior$U_hat[sel], table2_member = A_t2[sel], A_size = sum(A_t2))
  })
  summarise_oc(reps, mc$J, mc$K, true_pE, true_pT, mc$deltaE, mc$deltaT)
}

RUNNERS <- list(
  "Comb-BOIN12" = function(mc, tE, tT, seeds) run_oc_comb_boin12(mc, tE, tT, seeds),
  "EffTox"      = function(mc, tE, tT, seeds) run_oc_efftox(mc, tE, tT, seeds, fast = TRUE),
  "EffTox-TITE" = function(mc, tE, tT, seeds) run_oc_efftox_tite(mc, tE, tT, seeds, fast = TRUE),
  "COMIC"       = function(mc, tE, tT, seeds) run_oc_comic(mc, tE, tT, seeds),
  "uTPI-Comb"   = function(mc, tE, tT, seeds) run_oc_utpi_comb(mc, tE, tT, seeds)
)

sensitivity_note <- function(run_fun, mc, true_pE, true_pT, n_sens = min(4, mc$Nsim)) {
  mc_sens <- mc; mc_sens$wT <- mc$wT * 1.5
  oc_sens <- tryCatch(run_fun(mc_sens, true_pE, true_pT, 9000 + seq_len(n_sens)), error = function(e) NULL)
  if (is.null(oc_sens)) return("sensitivity re-run failed")
  sprintf("increasing wT by 50%% changed the modal cell to (%s,%s) and Pr(select unsafe) to %s [n=%d replicates]",
          oc_sens$modal_cell[1], oc_sens$modal_cell[2],
          ifelse(is.na(oc_sens$Pr_select_unsafe), "NA", sprintf("%.2f", oc_sens$Pr_select_unsafe)), n_sens)
}

report_design <- function(design_name, scenario_name, scenario, mc, oc, note) {
  cat("\n------------------------------------------------------------\n")
  cat("DESIGN:", design_name, " | SCENARIO:", scenario_name, "\n")
  cat("------------------------------------------------------------\n")
  cat("(J, K) = (", mc$J, ",", mc$K, ")\n")
  cat("(aE, bE) = (", mc$aE, ",", mc$bE, ")   (aT, bT) = (", mc$aT, ",", mc$bT, ")\n")
  cat("(deltaE, deltaT, piE, piT) = (", mc$deltaE, ",", mc$deltaT, ",", mc$piE, ",", mc$piT, ")\n")
  cat("(wE, wT) = (", mc$wE, ",", mc$wT, ")\n")
  cat("n =", mc$n, "   Nsim =", mc$Nsim, "\n")
  cat("Source of the truth surfaces:", scenario$source, "\n\n")
  
  if (all(!is.na(oc$modal_cell))) {
    j <- oc$modal_cell[1]; k <- oc$modal_cell[2]
    cat("Modal cell across replicates [most frequently part of the recommended set] = (", j, ",", k, ")\n")
    if (oc$n_tied_modal > 1) {
      cat("  (", oc$n_tied_modal, "cells tied on selection frequency; lowest-dose tie-break applied for this display cell)\n")
    }
    cat("  Ehat =", round(oc$Ehat_at_modal, 3), " That =", round(oc$That_at_modal, 3),
        " U =", round(oc$Uhat_at_modal, 3), "  (averaged across every recommended cell in every replicate)\n")
    cat("  Membership in A: satisfied for", round(100 * oc$prop_modal_table2_member, 1),
        "% of recommended cells across all replicates\n\n")
  } else {
    cat("Recommended (j, k): no dose pair selected in any replicate (",
        round(100 * oc$prop_no_selection, 1), "% of trials had no admissible dose)\n\n")
  }
  
  cat("Selection-probability heatmap (rows = Agent A level j, cols = Agent B level k):\n")
  print(round(oc$heatmap, 2))
  cat("No dose selected in", round(100 * oc$prop_no_selection, 1), "% of replicates\n\n")
  
  cat("Pr(select T > deltaT) =", round(oc$Pr_select_unsafe, 3), "\n")
  cat("Pr(select E < deltaE) =", round(oc$Pr_select_ineffective, 3), "\n\n")
  if (!is.na(oc$prop_replicates_with_tie) && oc$prop_replicates_with_tie > 0) {
    cat("Replicates whose OWN final selection was itself a tie (>1 admissible cell shared the\n",
        "  max posterior utility within that single trial -- ALL such cells were counted, no\n",
        "  tie-break applied):", round(100 * oc$prop_replicates_with_tie, 1), "%\n\n", sep = "")
  }
  cat("Mean true E at selection =", round(oc$mean_true_E_at_selection, 3), "\n")
  cat("Mean true T at selection =", round(oc$mean_true_T_at_selection, 3), "\n\n")
  cat("Mean |A| =", round(oc$mean_A_size, 2), "\n\n")
  cat("Sensitivity summary:", note, "\n\n")
  cat("Assumptions: cells are modelled independently via Beta-Binomial updates\n",
      "(no explicit borrowing across the grid; EffTox-TITE is the exception, using\n",
      "a time-to-event weighted likelihood for continuous accrual -- see efftox.R),\n",
      "and utility is linear in (Ehat,That); demonstration truth surfaces should be\n",
      "replaced by calibrated scenario matrices informed by biology, PK/PD, or prior\n",
      "evidence, for substantive use. Ties are handled with two different policies:\n",
      "the TRUE OBDC and each design's own FINAL recommended OBDC are REPORTED\n",
      "CONCLUSIONS, not operational necessities, so they are NEVER tie-broken --\n",
      "every admissible cell sharing the max utility is reported. The 'modal cell\n",
      "across replicates' shown above is a different thing (a display summary of\n",
      "selection frequency over many simulated trials), and it IS tie-broken (lowest\n",
      "total standardized dose j+k) purely so the report has one cell to highlight;\n",
      "how often that display tie happened, and how often a single trial's own final\n",
      "selection was itself untied, are both reported above.\n", sep = "")
}

## =============================================================================
## run_obdc_comparison() -- the function to use once you SOURCE this file
## =============================================================================
## Sourcing this file (source("execution.R")) does NOT automatically run the
## full multi-scenario report -- it just loads the scenarios, the design
## engines, and this function, so you can call it yourself with whichever
## scenario and arguments you want. (Running `Rscript execution.R` from a
## terminal still gives you the full default report, as before -- see the
## RUN section at the bottom.)
##
## scenario: either the NAME of one of the fixed scenarios (a character
##   string -- see names(ALL_SCENARIOS) for the full list, e.g. "Scenario 1"),
##   or a scenario you built yourself with make_custom_scenario(J, K, true_pT,
##   true_pE) from scenarios.R. Either way, J and K are taken directly from
##   the scenario -- you never need to type them in separately.
## ...: any of the Table 2 arguments (deltaE, deltaT, piE, piT, wE, wT, aE,
##   bE, aT, bT, n, Nsim) -- edit only the ones you want to change from the
##   defaults at the top of this file.
## Returns all five designs' operating-characteristics results (invisibly),
## and -- unless verbose = FALSE -- prints the same report as execution.R's
## default run, for all five designs on the one scenario you chose.
run_obdc_comparison <- function(scenario,
                                deltaE = DELTA_E, deltaT = DELTA_T,
                                piE = PI_E, piT = PI_T,
                                wE = W_E, wT = W_T,
                                aE = A_E, bE = B_E, aT = A_T, bT = B_T,
                                n = N_PER_CELL, Nsim = N_SIM,
                                verbose = TRUE) {
  
  ## Resolve the scenario -- a name (character) or a ready-made scenario list.
  if (is.character(scenario)) {
    if (!scenario %in% names(ALL_SCENARIOS)) {
      stop("Unknown scenario name: '", scenario, "'.\nValid names are: ",
           paste(names(ALL_SCENARIOS), collapse = ", "), call. = FALSE)
    }
    scenario_label <- scenario
    sc <- ALL_SCENARIOS[[scenario]]
  } else if (is.list(scenario) && all(c("J", "K", "true_pT", "true_pE") %in% names(scenario))) {
    scenario_label <- if (!is.null(scenario$source)) scenario$source else "Custom scenario"
    sc <- scenario
  } else {
    stop("scenario must be either the name of a fixed scenario (see names(ALL_SCENARIOS)) ",
         "or a scenario list with J, K, true_pT, true_pE (e.g. from make_custom_scenario()).",
         call. = FALSE)
  }
  
  ## The common Table-2 parameter set. J and K are extracted straight from
  ## the scenario -- the dose grid belongs to the scenario, not something you
  ## set independently.
  mc <- list(deltaE = deltaE, deltaT = deltaT, piE = piE, piT = piT,
             wE = wE, wT = wT, aE = aE, bE = bE, aT = aT, bT = bT,
             n = n, Nsim = Nsim, J = sc$J, K = sc$K)
  seeds <- seq_len(mc$Nsim)
  
  if (verbose) report_true_obdc(scenario_label, sc, mc)
  
  results <- list()
  for (design_name in names(RUNNERS)) {
    if (verbose) cat("Running", design_name, "on", scenario_label, "...\n")
    oc <- RUNNERS[[design_name]](mc, sc$true_pE, sc$true_pT, seeds)
    note <- sensitivity_note(RUNNERS[[design_name]], mc, sc$true_pE, sc$true_pT)
    if (verbose) report_design(design_name, scenario_label, sc, mc, oc, note)
    results[[design_name]] <- oc
  }
  invisible(results)
}

## -----------------------------------------------------------------------------
## SIMPLE EXAMPLE -- safe to copy/paste into your R console after sourcing
## this file (source("execution.R")):
##
##   # All 5 designs on the paper's "Scenario 1", using the default Table-2
##   # arguments. J and K are NOT typed in by hand -- they come straight from
##   # the scenario:
##   results <- run_obdc_comparison("Scenario 1")
##
##   # Same scenario, but edit whichever arguments you want:
##   results <- run_obdc_comparison("Scenario 1", deltaE = 0.25, wT = 1.5, Nsim = 20)
##
##   # Any other fixed scenario -- see names(ALL_SCENARIOS) for the full list:
##   results <- run_obdc_comparison("Scenario 9", piE = 0.75)
##
##   # A scenario you built yourself, instead of a name:
##   my_scn <- make_custom_scenario(J = 4, K = 3,
##                                   true_pT = matrix(c(0.05,0.10,0.15,0.20,
##                                                       0.10,0.15,0.20,0.25,
##                                                       0.15,0.20,0.25,0.30),
##                                                     nrow = 4, byrow = TRUE),
##                                   true_pE = matrix(c(0.10,0.20,0.30,0.35,
##                                                       0.20,0.30,0.40,0.45,
##                                                       0.15,0.25,0.35,0.40),
##                                                     nrow = 4, byrow = TRUE))
##   results <- run_obdc_comparison(my_scn)
## -----------------------------------------------------------------------------

## =============================================================================
## RUN  --  only fires automatically for `Rscript execution.R` from a
## terminal (non-interactive). If you `source("execution.R")` from an
## interactive R/RStudio session instead, this block is skipped -- you get
## run_obdc_comparison() and everything it needs, ready to call yourself; see
## the SIMPLE EXAMPLE above.
## =============================================================================

if (!interactive()) {
  
  missing_scn <- setdiff(SCENARIOS_TO_RUN, names(ALL_SCENARIOS))
  if (length(missing_scn) > 0) {
    stop("Unknown scenario name(s) in SCENARIOS_TO_RUN: ", paste(missing_scn, collapse = ", "),
         "\nValid names are: ", paste(names(ALL_SCENARIOS), collapse = ", "))
  }
  
  cat("=========================================================\n")
  cat(" OBDC design comparison  --  ", length(SCENARIOS_TO_RUN), "scenario(s) x", length(RUNNERS), "designs\n")
  cat("=========================================================\n")
  
  for (scn_name in SCENARIOS_TO_RUN) {
    scenario <- ALL_SCENARIOS[[scn_name]]
    mc <- CFG; mc$J <- scenario$J; mc$K <- scenario$K
    seeds <- seq_len(mc$Nsim)
    
    cat("\n\n=========================================================\n")
    cat("SCENARIO:", scn_name, " (source:", scenario$source, ")\n")
    cat("=========================================================\n")
    report_true_obdc(scn_name, scenario, mc)
    
    for (design_name in names(RUNNERS)) {
      cat("Running", design_name, "on", scn_name, "...\n")
      oc <- RUNNERS[[design_name]](mc, scenario$true_pE, scenario$true_pT, seeds)
      note <- sensitivity_note(RUNNERS[[design_name]], mc, scenario$true_pE, scenario$true_pT)
      report_design(design_name, scn_name, scenario, mc, oc, note)
    }
  }
  
  cat("\nDone. These are original re-implementations built from the source papers'\n",
      "equations, not the original authors' validated software -- treat as a\n",
      "design-exploration tool, not validated trial software.\n", sep = "")
  cat("\nFor an interactive point-and-click version of these same comparisons, see\n",
      "the hosted app linked in README.md.\n", sep = "")
  
}