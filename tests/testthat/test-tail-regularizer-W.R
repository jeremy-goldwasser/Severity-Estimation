tail_regularizer_mat_old <- function(n, delay_distr) {
  n_diffs <- n - 1
  d <- length(delay_distr) - 1
  delay_cdf <- cumsum(delay_distr) # P(die within k most recent days)
  reg_weights <- 1/sqrt(delay_cdf)
  W <- matrix(0, n_diffs, n_diffs)
  indices <- seq(n_diffs - d, n_diffs)
  W[cbind(indices, indices)] <- rev(reg_weights)
  W
}


test_that("tail_regularizer_mat matches tail_regularizer_mat_old", {
  # Example parameters
  n <- 6
  delay <- c(0.2, 0.3, 0.5)
  locs <- 1:n
  # Run both
  W_new <- tail_regularizer_mat(n, delay, locs)
  W_old <- tail_regularizer_mat_old(n, delay)
  
  expect_equal(as.matrix(W_new), W_old)
  
  n <- 8
  delay <- c(0.1, 0.2, 0.3, 0.4)
  locs <- 1:n
  
  W_new <- tail_regularizer_mat(n, delay, locs)
  W_old <- tail_regularizer_mat_old(n, delay)
  
  expect_equal(as.matrix(W_new), W_old)
})

test_that("tail_regularizer_mat handles missing delay locations correctly", {
  n <- 10
  delay <- c(0.1, 0.2, 0.3, 0.4)  # only first 4 days have delay mass
  locs <- c(1:4, 6, 8, 10)          # skip some days entirely
  
  # Expected behaviour:
  # 1. Pad delay with zeros to the full length for missing days.
  # 2. Keep only entries at `locs`.
  # 3. Normalize to sum to 1.
  # 3. Replace positive entries with 1/sqrt(cumsum(positive entries in order)).
  full_delay <- c(delay, rep(0, max(locs) - length(delay)))
  partial_delay <- rev(rev(full_delay)[locs])
  partial_delay <- partial_delay / sum(partial_delay)
  pos_vals <- partial_delay[partial_delay > 0]
  partial_delay[partial_delay > 0] <- 1 / sqrt(cumsum(pos_vals))
  
  expected_diag <- rev(partial_delay)[-1]
  
  W <- tail_regularizer_mat(n, delay, locs)
  diag_W <- diag(as.matrix(W))
  
  expect_equal(diag_W, expected_diag)
})
