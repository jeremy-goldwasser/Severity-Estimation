test_that("smoother refactor matches the old version", {
  x <- rnorm(100)
  zec <- smooth_counts_old(x, 6, TRUE)
  zer <- smooth_counts_old(x, 6, FALSE)
  zoc <- smooth_counts_old(x, 7, TRUE)
  zor <- smooth_counts_old(x, 7, FALSE)

  # Even window
  expect_equal(smooth_counts(x, 6, TRUE), zec)
  expect_equal(smooth_counts(x, 6, FALSE), zer)

  # Odd window
  expect_equal(smooth_counts(x, 7, TRUE), zoc)
  expect_equal(smooth_counts(x, 7, FALSE), zor)
})


test_that("seven day smoothing refactor matches the old version", {
  x <- rnorm(100)
  zc <- seven_day_smoothing_old(x, TRUE)
  zr <- seven_day_smoothing_old(x, FALSE)
  expect_equal(seven_day_smoothing(x, TRUE), zc)
  expect_equal(seven_day_smoothing(x, FALSE), zr)
})

test_that("seven day smoothing refactor matches the old version, named", {
  x <- rnorm(100)
  names(x) <- sample(letters[1:26], 100, TRUE)
  zc <- seven_day_smoothing_old(x, TRUE)
  zr <- seven_day_smoothing_old(x, FALSE)
  expect_equal(seven_day_smoothing(x, TRUE), zc)
  expect_equal(seven_day_smoothing(x, FALSE), zr)
})
