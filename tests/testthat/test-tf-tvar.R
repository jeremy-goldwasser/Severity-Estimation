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


delay <- make_delay_distr(10, 7, 60)
ndeath <- length(dat$dd)
tvar_delay <- outer(rep(1, length(dat$hh)), delay)


test_that("compute works (interface to Clarabel)", {
  Zs <- make_Z_mat_simple(ndeath, dat$hh, delay)
  Dmat <- make_tvar_delay_mat(tvar_delay, TRUE)
  Zt <- make_Z_mat_new(dat$hh, Dmat)
  delay_distr_T <- delay
  Y <- dat$dd
  maxmax <- max(max(abs(Zt)), max(abs(Zs)), abs(Y))

  cons <- compute_tf_hfrs(
    Zs / maxmax,
    dat$dd / maxmax,
    lambda = 100,
    order = 2,
    poisson = TRUE,
    real_time = FALSE,
    delay_distr_T = delay_distr_T
  )
  tvar <- compute_tf_hfrs(
    Zt / maxmax,
    dat$dd / maxmax,
    lambda = 100,
    order = 2,
    poisson = TRUE,
    real_time = FALSE,
    delay_distr_T = delay_distr_T
  )
  expect_equal(cons, tvar)

  cons <- compute_tf_hfrs(
    Zs / maxmax,
    dat$dd / maxmax,
    lambda = 100,
    order = 2,
    poisson = TRUE,
    real_time = TRUE,
    delay_distr_T = delay_distr_T
  )
  tvar <- compute_tf_hfrs(
    Zt / maxmax,
    dat$dd / maxmax,
    lambda = 100,
    order = 2,
    poisson = TRUE,
    real_time = TRUE,
    delay_distr_T = delay_distr_T
  )
  expect_equal(cons, tvar)
})

test_that("CV works", {
  cv_cons <- cv_lambda_tf(dat$hh, dat$dd, delay, lambdas = c(10, 100), K = 3)
  cv_tvar <- cv_lambda_tf(
    dat$hh,
    dat$dd,
    delay,
    tvar_delay = tvar_delay,
    lambdas = c(10, 100),
    K = 3
  )
  expect_equal(cv_cons$full_fit, cv_tvar$full_fit)
  expect_equal(cv_cons$cv_errors, cv_tvar$cv_errors)

  cv_cons <- cv_lambda_tf(
    dat$hh,
    dat$dd,
    delay,
    lambdas = c(10, 100),
    real_time = TRUE,
    K = 3
  )
  cv_tvar <- cv_lambda_tf(
    dat$hh,
    dat$dd,
    delay,
    tvar_delay = tvar_delay,
    lambdas = c(10, 100),
    real_time = TRUE,
    K = 3
  )
  expect_equal(cv_cons$full_fit, cv_tvar$full_fit)
  expect_equal(cv_cons$cv_errors, cv_tvar$cv_errors)
})

test_that("FV works", {
  fv_cons <- fv_gamma_tf(
    dat$hh,
    dat$dd,
    lambda = 10,
    gammas = c(10, 100),
    delay_distr = delay
  )
  fv_tvar <- fv_gamma_tf(
    dat$hh,
    dat$dd,
    lambda = 10,
    gammas = c(10, 100),
    tvar_delay = tvar_delay,
    delay_distr = delay
  )
  expect_equal(fv_cons, fv_tvar)
})

test_that("RT tuning works", {
  tuned_cons <- tune_trendfilter_rt(
    dat$hh,
    dat$dd,
    gammas = c(10, 100),
    delay_distr = delay,
    lambdas = 10,
    K = 3
  )
  tuned_tvar <- tune_trendfilter_rt(
    dat$hh,
    dat$dd,
    gammas = c(10, 100),
    delay_distr = delay,
    tvar_delay = tvar_delay,
    lambdas = 10,
    K = 3
  )
  expect_equal(tuned_cons$full_fit, tuned_tvar$full_fit)
  expect_equal(tuned_cons$fv_errors, tuned_tvar$fv_errors)
})
