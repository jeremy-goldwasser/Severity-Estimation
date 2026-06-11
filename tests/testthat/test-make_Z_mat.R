test_that("delay mat constructor produces the correct structure", {
  d1 <- make_delay_mat(6, 10, 1:5)
  
  hosps <- rep(1, 10)
  names(hosps) <- seq(as.Date("2021-01-01"), length.out = 10, by = "1 day")
  est_death_dates <- names(hosps)[-c(1:4)]
  Zold <- make_Z_mat_old(hosps, 1:5, as.Date(est_death_dates)) |>
    drop0()
  expect_identical(d1, Zold)
})


test_that("new Z mat constructor matches the old version (+ sparsity)", {
  d1 <- make_delay_mat(6, 10, 1:5 / 15)
  
  hosps <- rep(1, 10)
  names(hosps) <- seq(as.Date("2021-01-01"), length.out = 10, by = "1 day")
  est_death_dates <- names(hosps)[-c(1:4)]
  Zold <- make_Z_mat_old(hosps, 1:5 / 15, as.Date(est_death_dates)) |>
    drop0()
  Znew <- make_Z_mat(unname(hosps), d1)
  expect_identical(Zold, Znew)
  
  hosps <- c(2, 10, 3, 15, 6, 5, 3, 8, 7, 10)
  names(hosps) <- seq(as.Date("2021-01-01"), length.out = 10, by = "1 day")
  est_death_dates <- names(hosps)[-c(1:4)]
  Zold <- make_Z_mat_old(hosps, 1:5 / 15, as.Date(est_death_dates)) |>
    drop0()
  Znew <- make_Z_mat(unname(hosps), d1)
  expect_identical(Zold, Znew)
})

test_that("new Z mat constructor handles indices, renormalizes", {
  d1 <- make_delay_mat(6, 10, 1:5 / 15)
  
  hosps <- rep(1, 10) # removes the impact of the multiplication
  idc <- c(1:4, 6, 8, 10)
  Znew <- make_Z_mat(hosps, d1, idc = idc)
  ex <- d1[,idc]
  ex <- ex / rowSums(ex)
  expect_equal(Znew, ex)
  
  idr <- c(1, 3, 6)
  Znew <- make_Z_mat(hosps, d1, idr)
  ex <- d1
  expect_equal(Znew, ex[idr,])
  Znew <- make_Z_mat(hosps, d1, idr, idc)
  ex <- d1[,idc]
  ex <- ex / rowSums(ex)
  expect_equal(Znew, ex[idr,])
  
  hosps <- c(2, 10, 3, 15, 6, 5, 3, 8, 7, 10)
  Znew <- make_Z_mat(hosps, d1, idc = idc)
  ex1 <- ex %*% Diagonal(x = hosps[idc])
  expect_equal(Znew, ex1)
  
  Znew <- make_Z_mat(hosps, d1, idc = idc, idr = idr)
  expect_equal(Znew, ex1[idr,])
})

