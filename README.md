# OBDC-Comparison

This repository contains the simulation code used to conduct the analyses
for the following manuscript on Bayesian designs for identifying the
**optimal biological dose combination (OBDC)** — the dose pair achieving
sufficient biological activity with acceptable tolerability — in two-agent
Phase I/II oncology trials, in line with the FDA's Project Optimus
initiative. If you use this repository, please cite:

> Mukherjee A¹, Takeda K², Wason J M S¹. *Optimal Dose Combination Selection
> in Oncology Trials Using Bayesian Designs: Statistical, Regulatory, and
> Operational Insights* (manuscript in preparation).
>
> ¹Population Health Sciences Institute, Newcastle University, Newcastle upon
> Tyne, UK. ²Quantitative Science and Evidence Generation, Astellas Pharma
> Global Development Inc., Northbrook, USA.
>
> Corresponding author: Ayon Mukherjee ([ayon.mukherjee@newcastle.ac.uk](mailto:ayon.mukherjee@newcastle.ac.uk)).

It provides working R implementations of four representative designs spanning
the three principal statistical paradigms for OBDC estimation: a model-based
design (**EffTox**), a model-assisted design (**Comb-BOIN12**), and two
utility-integrated designs (**COMIC** and **uTPI-Comb**). Each design is
built from a common set of configurable parameters and can be run either as a
command-line report or through an interactive Shiny application, allowing
users to explore scenario assumptions and compare operating characteristics
across designs without needing to read the underlying statistical code.

This reporsitory contains Six self-contained .R files (the scenarios, the four design
engines, one "just-run-it" script and an R Shiny app makes up the bundle.


## Try the interactive app online

OBDC Compare App:
https://22c7ba-damitri-kundu.shinyapps.io/OBDC-Compare-App/
![Version](https://img.shields.io/badge/version-1.0-blue)

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

## References
1. Liang H, Yang Y N, Yuan M. uTPI-Comb: an optimal Bayesian dose-allocation method in two-agent phase I/II clinical trials. JUSTC, 2024, 54(12): 1206. DOI: 10.52396/JUSTC-2024-0104.

2. Lu M, Zhang J, Yuan Y, Lin R. Comb-BOIN12: A Utility-Based Bayesian Optimal Interval Design for Dose Optimization in Cancer Drug-Combination Trials. Statistics in Biopharmaceutical Research, 2025, 17(2): 266-276.

3. Thall P F, Cook J D. Dose-finding based on efficacy-toxicity trade-offs. Biometrics, 2004, 60(3): 684-693.
Chen K, Takeda K, Yuan Y. COMIC: A Bayesian dose optimization design for drug combination in multiple indications with application to CAR-T therapies. Statistics in Medicine, 2025, 44(10-12): e70107.

4. Mukherjee A¹, Takeda K², Wason J M S¹. Optimal Dose Combination Selection in Oncology Trials Using Bayesian Designs: Statistical, Regulatory, and Operational Insights (manuscript in preparation). ¹Population Health Sciences Institute, Newcastle University, Newcastle upon Tyne, UK. ²Quantitative Science and Evidence Generation, Astellas Pharma Global Development Inc., Northbrook, USA. Corresponding author: Ayon Mukherjee (ayon.mukherjee@newcastle.ac.uk).
