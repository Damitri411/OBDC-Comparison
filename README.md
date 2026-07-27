# OBDC Designs -- 7-file bundle

Seven self-contained files: the scenarios, the four design engines, one
"just-run-it" script, and one single-file Shiny app.

## Try the interactive app online

OBDC Compare App:
https://22c7ba-damitri-kundu.shinyapps.io/OBDC-Compare-App/

| File | What it is |
|------|------------|
| `scenarios.R` | Scenarios 1, 3, 5, 9, 10 from Table 1 of the uTPI-Comb paper (Liang, Yang & Yuan, 2024), plus `make_custom_scenario()` for building your own J x K grid. |
| `comb_boin12.R` | Comb-BOIN12 (model-assisted). Lu, Zhang, Yuan & Lin (2025). |
| `efftox.R` | EffTox (fully model-based). Thall & Cook (2004); Brock et al. (2017). |
| `comic.R` | COMIC (utility-integrated, multi-indication). Chen, Takeda & Yuan (2025). |
| `utpi_comb.R` | uTPI-Comb (utility-integrated, zone-based). Liang, Yang & Yuan (2024). |
| `execution.R` | **Run this for a plain-text report.** Edit the settings block at the top, then `Rscript execution.R`. |

## Quick start

**Command-line report** (no packages beyond base R needed):
```bash
Rscript execution.R
```
Everything you'd normally want to change — which scenarios to run, deltaE,
deltaT, piE, piT, wE, wT, the priors, n, Nsim, and an optional custom scenario
— is in one clearly marked block at the top of the file. Takes 2-3 minutes for
the default (all 5 scenarios x all 4 designs, Nsim=100).

Was tested end-to-end before delivery (execution.R run to completion).

**Or, from an R session** (`source("execution.R")`): this loads everything
without running the full report, and gives you `run_obdc_comparison()` to
call yourself, e.g.:
```r
source("execution.R")
results <- run_obdc_comparison(
  scenario = "Scenario 1",  # or any name in names(ALL_SCENARIOS), or your own via make_custom_scenario()
  deltaE   = 0.25,          # minimal efficacy threshold
  deltaT   = 0.30,          # maximal toxicity threshold
  piE      = 0.80,          # posterior requirement, efficacy
  piT      = 0.80,          # posterior requirement, toxicity
  wE       = 1,             # utility weight on efficacy
  wT       = 1.5,           # utility weight on toxicity
  aE       = 1, bE = 1,     # efficacy Beta prior (aE, bE)
  aT       = 1, bT = 1,     # toxicity Beta prior (aT, bT)
  n        = 6,             # OC per-cell sample size
  Nsim     = 100,           # OC replicates
  verbose  = TRUE           # print the full report as it runs
)
```
See the "SIMPLE EXAMPLE" comment block in execution.R for more.

#ted trial software — see the
   closing notes in its output for the fuller list of simplifications made
   in each engine.
## Caution 
In the R shiny App (https://22c7ba-damitri-kundu.shinyapps.io/OBDC-Compare-App/), Nsim higher than 500 can cause delays in execution.

