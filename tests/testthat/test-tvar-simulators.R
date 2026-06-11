source("Code/helper_functions.R")

delay <- make_delay_distr(3, 2.2, 10)
hfrs <- runif(31)
hosps <- rpois(31, 100)
dd <- seq(as.Date("2025-03-01"), length.out = 31, by = "1 day")
names(hosps) <- dd
names(hfrs) <- dd
death_dates <- tail(dd, 21)

test_that("poisson binomial matches the old version", {
  set.seed(12345)
  oldpb <- sim_pb_deaths_old(death_dates, hosps, hfrs, delay)
  set.seed(12345)
  newpb <- sim_pb_deaths(death_dates, hosps, hfrs, delay)

  expect_equal(newpb, oldpb)
})

test_that("beta binomial matches the old version", {
  set.seed(12345)
  oldbb <- sim_bb_deaths_old(death_dates, hosps, hfrs, delay, dispersion = 4)
  set.seed(12345)
  newbb <- sim_bb_deaths(death_dates, hosps, hfrs, delay, dispersion = 4)

  expect_equal(newpb, oldpb)
})
