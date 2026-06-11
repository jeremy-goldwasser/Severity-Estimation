# Severity-Estimation

Code and committed data for estimating time-varying disease severity — the
hospitalization-fatality ratio (HFR) — from noisy, right-censored surveillance
streams of hospitalizations and deaths.

The core method poses severity estimation as a Poisson trend-filtering
deconvolution: it jointly smooths the underlying signals and accounts for the
reporting delay between hospitalization and death, and is benchmarked against
simpler lagged-ratio estimators. Methods are evaluated on simulated problems
(retrospective, real-time, and time-varying-delay regimes) and on real U.S.
COVID-19 data pulled from the [Delphi Epidata API](https://cmu-delphi.github.io/delphi-epidata/).

This repository is figures-from-committed-data: every paper figure can be
regenerated from the `Data/batch-results/*.rds` files already checked in,
without re-running the compute-heavy experiments.

## Repository layout

| Path        | Contents                                                        |
| ----------- | --------------------------------------------------------------- |
| `Code/`     | Data download, core estimators, experiments, and analysis. See `Code/README.md`. |
| `Data/`     | `raw-data/` inputs, `Variants/` variant proportions, and the committed `batch-results/` used to build figures. |
| `Figures/`  | Generated figures, in `Color/` and `BW/` variants.              |
| `tests/`    | Unit tests for the core estimator functions.                   |

## Installation

Analyses are written in R. The main dependencies (all on CRAN unless noted):

- **Estimation / optimization:** `CVXR`, `clarabel`, `dspline`, `Matrix`, `VGAM`, `extraDistr`
- **Data access & wrangling:** `epidatr`, `epiprocess`, `tidyverse`, `dplyr`, `purrr`, `lubridate`, `zoo`, `RcppRoll`
- **Experiments & plotting:** `batchtools`, `ggplot2`, `here`, `cli`, `rlang`

The latest version of each should work. We recommend `conda` (e.g. the
`r-base` / `r-essentials` channels) to manage any dependency conflicts.

## Estimating severity from two incidence series

If you have a **primary** incidence series (e.g. hospitalizations) and a lagged
**secondary** series (e.g. deaths), the core method deconvolves the
reporting delay between them to estimate a time-varying severity rate:

```r
source("Code/helper_functions.R")
source("Code/estimate.R")

delay_distr <- make_delay_distr(Mean = 19, Sd = 17, d = 40)   # delay over lags 0..40

# Date-named vectors; primary leads secondary by d: length(hosps) == length(deaths) + 40
fit <- cv_lambda_tf(hosps, deaths, delay_distr, order = 2)
severity <- fit$full_fit[, fit$cvstats$i.min]                  # CV-optimal severity series
```

See [`Code/README.md`](Code/README.md) for inputs, the real-time and
time-varying-delay variants, and the ratio-method baseline.

## Reproducing the figures

To rebuild the paper figures from the committed results, run from the repository
root:

```r
source("Code/Analysis/manuscript-figs-color.R")   # writes Figures/Color/
source("Code/Analysis/manuscript-figs-bw.R")       # writes Figures/BW/
```

These read only from `Data/`, so no experiments need to be re-run. The full
pipeline — from raw data download through the experiments to the committed
results — is documented in [`Code/README.md`](Code/README.md).
