## =============================================================================
## scenarios.R
## Five of the ten scenarios from Table 1 of:
##   Liang H, Yang Y N, Yuan M. uTPI-Comb: an optimal Bayesian dose-allocation
##   method in two-agent phase I/II clinical trials. JUSTC, 2024, 54(12): 1206.
##   DOI: 10.52396/JUSTC-2024-0104
##
## Table 1 reports each scenario as a 3 (Drug B levels) x 5 (Drug A levels)
## grid of TRUE toxicity and efficacy probabilities. We transcribe each matrix
## exactly as printed (rows = Drug B level 1..3, columns = Drug A level 1..5)
## and then transpose so that, consistent with every other file in this
## project, true_pT/true_pE are J x K matrices with rows = Agent A dose level
## (J = 5) and columns = Agent B dose level (K = 3).
##
## The paper defines: OBDC = the dose combination with the highest efficacy
## probability among those "deemed safe"; TDC = a dose combination among the
## safe doses with efficacy probability >= 45%. We reproduce the OBDC/TDC
## flags exactly as bolded/underlined in Table 1, for reference display only
## -- any "true admissible set" / "true OBDC" computed elsewhere in this
## project uses the USER'S chosen (deltaE, deltaT), rather than the paper's
## own (phiT = 0.35, phiE = 0.2) simulation settings, so it will not always
## coincide with the paper's bolded cells unless deltaT=0.35, deltaE=0.2.
##
## Also provides make_custom_scenario(), a small helper for building an
## arbitrary J x K scenario (used by execution.R and app.R for user-defined
## grids) in the same list structure as the five fixed scenarios below.
## =============================================================================

.mk <- function(rows_by_dosebrace) {
  # helper: build a 3x5 matrix (Drug B rows x Drug A cols) then transpose to 5x3
  m <- matrix(rows_by_dosebrace, nrow = 3, byrow = TRUE)
  t(m)
}

UTPI_PAPER_SCENARIOS <- list()

UTPI_PAPER_SCENARIOS[["Scenario 1"]] <- list(
  J = 5, K = 3,
  source = "scenario assumptions",
  paper_note = "TDC coincides with OBDC in this scenario (Table 1).",
  true_pT = .mk(c(0.05, 0.15, 0.30, 0.45, 0.55,
                  0.15, 0.35, 0.45, 0.55, 0.65,
                  0.35, 0.45, 0.55, 0.65, 0.75)),
  true_pE = .mk(c(0.05, 0.25, 0.50, 0.55, 0.60,
                  0.25, 0.50, 0.55, 0.60, 0.65,
                  0.50, 0.55, 0.60, 0.65, 0.70)),
  paper_OBDC = c(j = 3, k = 1)   # Drug A=3, Drug B=1 (0.50 efficacy, bolded+underlined cell)
)

UTPI_PAPER_SCENARIOS[["Scenario 3"]] <- list(
  J = 5, K = 3,
  source = "scenario assumptions",
  paper_note = "More TDC than OBDC dose combinations (Table 1).",
  true_pT = .mk(c(0.05, 0.10, 0.18, 0.25, 0.42,
                  0.10, 0.15, 0.23, 0.42, 0.43,
                  0.15, 0.23, 0.45, 0.50, 0.55)),
  true_pE = .mk(c(0.30, 0.45, 0.60, 0.45, 0.26,
                  0.20, 0.28, 0.45, 0.26, 0.18,
                  0.10, 0.14, 0.24, 0.18, 0.10)),
  paper_OBDC = c(j = 3, k = 1)   # 0.60 efficacy cell
)

UTPI_PAPER_SCENARIOS[["Scenario 5"]] <- list(
  J = 5, K = 3,
  source = "scenario assumptions",
  paper_note = "More TDC than OBDC dose combinations (Table 1).",
  true_pT = .mk(c(0.15, 0.21, 0.30, 0.42, 0.44,
                  0.24, 0.30, 0.42, 0.44, 0.51,
                  0.30, 0.33, 0.44, 0.51, 0.55)),
  true_pE = .mk(c(0.20, 0.45, 0.33, 0.15, 0.05,
                  0.35, 0.60, 0.45, 0.20, 0.15,
                  0.20, 0.45, 0.30, 0.15, 0.10)),
  paper_OBDC = c(j = 2, k = 2)   # 0.60 efficacy cell (Drug A=2, Drug B=2)
)

UTPI_PAPER_SCENARIOS[["Scenario 9"]] <- list(
  J = 5, K = 3,
  source = "scenario assumptions",
  paper_note = "More TDC than OBDC dose combinations; efficacy peak at LOW Agent-A dose (Table 1).",
  true_pT = .mk(c(0.15, 0.21, 0.30, 0.42, 0.44,
                  0.24, 0.30, 0.42, 0.44, 0.51,
                  0.30, 0.33, 0.44, 0.51, 0.55)),
  true_pE = .mk(c(0.21, 0.15, 0.12, 0.09, 0.05,
                  0.60, 0.45, 0.33, 0.24, 0.21,
                  0.45, 0.31, 0.24, 0.21, 0.17)),
  paper_OBDC = c(j = 1, k = 2)   # 0.60 efficacy cell (Drug A=1, Drug B=2)
)

UTPI_PAPER_SCENARIOS[["Scenario 10"]] <- list(
  J = 5, K = 3,
  source = "scenario assumptions",
  paper_note = "No TDC in this scenario -- toxicity is high everywhere, illustrating a case where the paper's own OBDC definition differs from the other scenarios (Table 1).",
  true_pT = .mk(c(0.50, 0.56, 0.65, 0.68, 0.72,
                  0.55, 0.62, 0.70, 0.72, 0.80,
                  0.60, 0.67, 0.75, 0.79, 0.85)),
  true_pE = .mk(c(0.52, 0.62, 0.70, 0.76, 0.79,
                  0.55, 0.66, 0.74, 0.79, 0.82,
                  0.58, 0.70, 0.78, 0.82, 0.85)),
  paper_OBDC = c(j = 1, k = 1)   # lowest dose, since toxicity is high everywhere
)

rm(.mk)

for (nm in names(UTPI_PAPER_SCENARIOS)) {
  sc <- UTPI_PAPER_SCENARIOS[[nm]]
  stopifnot(nrow(sc$true_pT) == sc$J, ncol(sc$true_pT) == sc$K,
            nrow(sc$true_pE) == sc$J, ncol(sc$true_pE) == sc$K)
}

## -----------------------------------------------------------------------------
## Helper: build a user-defined scenario in the same structure as the fixed
## ones above (J,K,true_pT,true_pE,source). J,K>=2. If true_pT/true_pE are not
## supplied, fills a smooth toy default (monotone toxicity, plateauing
## efficacy) purely as an editable starting point.
## -----------------------------------------------------------------------------
make_custom_scenario <- function(J, K, true_pT = NULL, true_pE = NULL,
                                 source = "scenario assumptions (user-defined)") {
  stopifnot(J >= 2, K >= 2)
  if (is.null(true_pT)) {
    true_pT <- outer(1:J, 1:K, function(a, b) 0.05 + 0.08 * (a - 1) + 0.08 * (b - 1))
    true_pT <- round(pmin(true_pT, 0.9), 2)
  }
  if (is.null(true_pE)) {
    true_pE <- outer(1:J, 1:K, function(a, b) 0.10 + 0.10 * (a + b - 2))
    true_pE <- round(pmin(true_pE, 0.75), 2)
  }
  stopifnot(nrow(true_pT) == J, ncol(true_pT) == K, nrow(true_pE) == J, ncol(true_pE) == K)
  list(J = J, K = K, true_pT = true_pT, true_pE = true_pE, source = source)
}

