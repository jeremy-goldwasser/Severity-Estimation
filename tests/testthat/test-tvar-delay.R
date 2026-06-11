library(tidyverse)
source("Code/batch-funs.R")
source("Code/funs.R")

raw_data <- read_rds(here::here("Data", "raw-data", "finalized.rds")) |>
  select(-deaths)
gt_hfrs <- read_rds(here::here("Data", "raw-data", "daily-hfrs.rds"))
meta_data <- read_rds(here::here("Data", "raw-data", "all-meta-data.rds"))

regions <- sort(meta_data$geo_value) # Includes US

## Global parameters
d <- 60L
ws <- c(1, 7, 15, 21, 29)
ord <- 0:2

## Delay mean modifier for misspecified retrospective case.
delay_mean_modifier <- c(0, 0, 0, -3, -2, -1, 1, 2, 3)
delay_sds <- c(.9, 1, .8, rep(.9, 6))

data_list <- enlist(
  raw_data,
  gt_hfrs,
  meta_data,
  delay_mean_modifier,
  delay_sds,
  ws,
  d
)

deaths <- create_regional_problem_retro(data_list, "ca") |>
  rename(deaths = sim_deaths)
dat <- data_list$raw_data |>
  filter(geo_value == "ca") |>
  left_join(deaths, by = join_by(geo_value, time_value)) |>
  process_tibble_to_named_vctrs()


test_that("tvar delay constructor works", {
  Pi <- outer(1:10, 4:0 / 5, "+")
  expect_equal(
    make_tvar_delay_mat(Pi, FALSE),
    t(bandSparse(10, 10, k = 0:4, diagonals = Pi))
  )
  expect_equal(
    spread_tvar_delay(Pi),
    bandSparse(10, 14, 0:4, diagonals = Pi)[-(1:4), -(1:4)]
  )
})

test_that("retro conv ratio is same with `delay_distr` or `tvar_delay`", {
  delay_distn <- make_delay_distr(10, 5, 60) #0:60
  tvar_delay <- outer(rep(1, length(dat$hh)), delay_distn)
  dvec <- compute_conv_hfrs_retro(dat$hh, dat$dd, delay_distn)
  tvec <- compute_conv_hfrs_retro(dat$hh, dat$dd, delay_distn, tvar_delay)
  expect_equal(dvec, tvec)
})


test_that("realtime conv ratio is same with `delay_distr` or `tvar_delay`", {
  delay_distn <- make_delay_distr(10, 5, 60) #0:60
  tvar_delay <- outer(rep(1, length(dat$hh)), delay_distn)
  dvec <- compute_conv_hfrs_realtime(
    dat$hh,
    dat$dd,
    delay_distn,
    dat$death_dates
  )
  tvec <- compute_conv_hfrs_realtime(
    dat$hh,
    dat$dd,
    delay_distn,
    dat$death_dates,
    tvar_delay
  )
  expect_equal(dvec, tvec)
})

test_that("realtime lagged ratio is same with `delay_distr` or `tvar_delay`", {
  lag <- 7
  tvar_lag <- rep(7, length(dat$hh))

  dvec <- compute_lagged_hfrs_fast(dat$hh, dat$dd, lag)
  tvec <- compute_lagged_hfrs_fast(dat$hh, dat$dd, -7, tvar_lag = tvar_lag)
  expect_equal(dvec, tvec)

  delay_distn <- make_delay_distr(10, 5, 60) #0:60
  tvar_delay <- outer(rep(1, length(dat$hh)), delay_distn)
  tlagged <- tune_lagged_hfrs_rt(
    dat$hh,
    dat$dd,
    as.Date(names(tail(dat$dd, 1))),
    delay_distn,
    c(1, 7, 15, 21, 29, 35),
    lag = lag
  )
  dlagged <- tune_lagged_hfrs_rt(
    dat$hh,
    dat$dd,
    as.Date(names(tail(dat$dd, 1))),
    delay_distn,
    c(1, 7, 15, 21, 29, 35),
    lag = lag,
    tvar_lag = tvar_lag,
    tvar_delay = tvar_delay
  )
  expect_equal(tlagged$full_fit, dlagged$full_fit)
  expect_equal(tlagged$fv_errors, dlagged$fv_errors)
})


test_that("retro lagged ratio is same with `delay_distr` or `tvar_delay`", {
  lag <- 7
  tvar_lag <- rep(7, length(dat$hh))
  dvec <- compute_lagged_hfrs_fast(dat$hh, dat$dd, lag, real_time = FALSE)
  tvec <- compute_lagged_hfrs_fast(
    dat$hh,
    dat$dd,
    -7,
    tvar_lag = tvar_lag,
    real_time = FALSE
  )
  expect_equal(dvec, tvec)
  delay_distn <- make_delay_distr(10, 5, 60) #0:60
  tvar_delay <- outer(rep(1, length(dat$hh)), delay_distn)
  tlagged <- tune_lagged_hfrs_retro(
    dat$hh,
    dat$dd,
    delay_distn,
    c(1, 7, 15, 21, 29, 35),
    lag = lag
  )
  dlagged <- tune_lagged_hfrs_retro(
    dat$hh,
    dat$dd,
    delay_distn,
    c(1, 7, 15, 21, 29, 35),
    lag = lag,
    tvar_lag = tvar_lag,
    tvar_delay = tvar_delay
  )
  expect_equal(tlagged$full_fit, dlagged$full_fit)
  expect_equal(tlagged$fv_errors, dlagged$fv_errors)
})

test_that("lagged ratio computes something reasonable with time varying lags", {
  n <- length(dat$hh)
  tvar_lag <- c(rep(7, n / 4), rep(5, n / 4), rep(9, n / 4), rep(7, n / 4))
  expect_no_condition(
    compute_lagged_hfrs_fast(
      dat$hh,
      dat$dd,
      -7,
      w = 21,
      tvar_lag = tvar_lag,
      real_time = FALSE
    )
  )
})
