test_that("Faster conv retro estimation matches the old version", {
  n <- 141
  hosps <- rpois(n, 40)
  names(hosps) <- seq(as.Date("2020-01-01"), length.out = n, by = "1 day")
  deaths <- rpois(n - 60, 20)
  names(deaths) <- names(hosps)[-c(1:60)]
  delay_distn <- make_delay_distr(12, 10, 60)
  old <- compute_conv_hfrs_retro_old(hosps, deaths, delay_distn, as.Date(names(deaths)[1:20]))
  new <- compute_conv_hfrs_retro(hosps, deaths, delay_distn)
  expect_equal(old, new)
})
