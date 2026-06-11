# Code

This directory holds the full analysis: data download, the core estimators, the
experiments that evaluate them, and the scripts that turn results into figures.

## Pipeline

The pipeline runs raw data → results → figures. Because the result files in
`Data/batch-results/` are committed, you can jump straight to step 4 to
reproduce the figures; steps 1–3 regenerate those results from scratch.

**1. Download data** — `download-data.R` pulls hospitalization, death, and
variant data from the Delphi Epidata API and writes the `Data/raw-data/*.rds`
inputs.

**2. Run experiments** — the scripts in `Experiments/` run trend filtering and
the ratio methods over the simulated and real problems via `batchtools`,
creating job registries. There is one script per evaluation regime:

| Script               | Regime                                                     |
| -------------------- | ---------------------------------------------------------- |
| `retro_batch.R`      | Retrospective (fully-observed) estimation on simulations.  |
| `rt_batch.R`         | Real-time estimation under right-censored reporting.       |
| `tvar_batch.R`       | Simulations with a time-varying reporting delay.           |
| `real-data_batch.R`  | Real U.S. COVID-19 hospitalization/death data.             |

**3. Collect & post-process** — the scripts in `Analysis/` gather and
post-process those registries into the committed `Data/batch-results/*.rds`:

- `collect-{retro,rt,time-var,real-data}.R` — collect each regime's registry.
- `process-real-data-results.R` — post-process the real-data results.
- `examine-retro-timing.R` — runtime comparisons.
- `examine-tail-regularization.R` — tail-regularization comparison.
- `maes-time-var.R` — time-varying-delay MAE summary.
- `collection-funs.R` — shared helpers for the above.

**4. Build figures** — `Analysis/manuscript-figs-color.R` and
`Analysis/manuscript-figs-bw.R` build every paper figure from the committed
results, writing to `Figures/Color/` and `Figures/BW/` respectively. These read
only from `Data/`, so the experiments do not need to be re-run.

## Core method / library code

- `estimate.R` — Poisson trend-filtering deconvolution and the ratio methods.
- `batch-funs.R` — per-region problem setup and the estimator drivers.
- `funs.R`, `helper_functions.R` — data processing, simulation, and utilities.

## Using the estimator

The method estimates a time-varying severity rate (a secondary-to-primary
ratio) from two incidence series — a **primary** signal (e.g. hospitalizations)
and a lagged **secondary** signal (e.g. deaths) — given the reporting-delay
distribution between them. The estimator deconvolves the delay rather than
applying a fixed lag, so it recovers severity even when the delay is long or the
signals are changing quickly.

Inputs:

- `hosps`, `deaths` — date-named numeric vectors of the primary and secondary
  incidence. The primary series must **lead** the secondary by `d` days, where
  `d` is the maximum modeled delay: `length(hosps) == length(deaths) + d`.
- `delay_distr` — the delay distribution over lags `0..d`, built with
  `make_delay_distr(Mean, Sd, d)` (a discretized gamma).

```r
source("Code/helper_functions.R")   # make_delay_distr, smoothing, utilities
source("Code/estimate.R")           # the estimators (loads required libraries)

# Delay from a primary event to the secondary event: mean 19 days, sd 17,
# truncated at d = 40 days.
d <- 40
delay_distr <- make_delay_distr(Mean = 19, Sd = 17, d = d)

# Primary leads secondary by d:  length(hosps) == length(deaths) + d
hosps  <- set_names(primary_counts,   as.character(primary_dates))
deaths <- set_names(secondary_counts, as.character(secondary_dates))

# --- Trend-filtering deconvolution (the core method) ---
fit <- cv_lambda_tf(hosps, deaths, delay_distr, order = 2)   # cross-validates lambda
severity <- fit$full_fit[, fit$cvstats$i.min]                # series at the CV-optimal lambda
# fit$full_fit has one column per candidate lambda; cvstats$i.min picks the
# CV minimizer (cvstats$i.1se gives the 1-standard-error choice).

# --- Convolution-ratio baseline (no tuning) ---
severity_ratio <- compute_conv_hfrs_retro(hosps, deaths, delay_distr)
```

Both return a severity series named by the secondary-signal dates it can
support (the first `d` dates are consumed by the delay window).

Variants of the call:

- **Real-time / right-censored** estimation: pass `real_time = TRUE` (and a
  `gamma > 0` tail penalty) to `cv_lambda_tf` so recent, still-incomplete
  secondary counts are handled — see `tune_trendfilter_rt` and `rt_batch.R`.
- **Time-varying delay**: pass a `tvar_delay` matrix (one delay row per primary
  date) instead of a single `delay_distr` — see `tvar_batch.R`.
- **Ratio method with window tuning**: `tune_conv_hfrs_retro(hosps, deaths,
  delay_distr, ws)` cross-validates the smoothing window `ws`.

For complete, runnable end-to-end examples, see the regime drivers in
`Experiments/` and the unit tests in `tests/testthat/`.

## Tests

`tests/` exercises the core functions in `estimate.R`, `batch-funs.R`,
`funs.R`, and `helper_functions.R`.
