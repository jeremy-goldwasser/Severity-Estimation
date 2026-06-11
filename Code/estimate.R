library(CVXR)
library(lubridate)
library(Matrix)
library(dplyr)
library(purrr)
library(zoo)
library(dspline)
library(rlang)
library(RcppRoll)
library(cli)


# Lagged HFR ------------------------------------------------------------------
# compute_lagged_hfr <- function(
#   hosps,
#   deaths,
#   t,
#   lag = 19,
#   w = 1,
#   real_time = TRUE
# ) {
#   w <- ifelse(w == 0, 1, w)
#   if (is.Date(date) == FALSE) {
#     t <- as.Date(t)
#   }

#   centered <- ifelse(real_time, FALSE, TRUE)
#   n_into_past <- ifelse(centered, floor((w - 1) / 2), w - 1)
#   n_into_future <- ifelse(centered, max(0, floor(w / 2)), 0)
#   if (real_time == FALSE) {
#     # Y_{t+l}/X_{t}
#     hosps_around_t <- hosps[
#       names(hosps) >= t - n_into_past & names(hosps) <= t + n_into_future
#     ]
#     deaths_in_around_l <- deaths[
#       names(deaths) >= t + lag - n_into_past &
#         names(deaths) <= t + lag + n_into_future
#     ]
#     if (length(hosps_around_t) == 0 | length(deaths_in_around_l) == 0) {
#       return(NA)
#     }
#     hosps_t <- sum(hosps_around_t)
#     deaths_in_l <- sum(deaths_in_around_l)
#     return(deaths_in_l / hosps_t)
#   }
#   # Y_{t}/X_{t-l}; different strategies for window
#   hosps_around_l_ago <- hosps[
#     names(hosps) >= t - n_into_past - lag &
#       names(hosps) <= t + n_into_future - lag
#   ]
#   deaths_around_t <- deaths[
#     names(deaths) >= t - n_into_past & names(deaths) <= t + n_into_future
#   ]
#   if (length(hosps_around_l_ago) == 0 | length(deaths_around_t) == 0) {
#     return(NA)
#   }
#   hosps_l_ago <- sum(hosps_around_l_ago)
#   deaths_t <- sum(deaths_around_t)
#   return(deaths_t / hosps_l_ago)
# }

# compute_lagged_hfrs <- function(
#   hosps,
#   deaths,
#   lag,
#   w = 1,
#   real_time = TRUE,
#   dates = NULL
# ) {
#   if (is.null(dates)) {
#     dates <- names(deaths)
#   }
#   lagged_hfrs <- sapply(
#     dates,
#     compute_lagged_hfr,
#     lag = lag,
#     w = w,
#     real_time = real_time,
#     hosps = hosps,
#     deaths = deaths
#   )
#   names(lagged_hfrs) <- dates
#   return(lagged_hfrs)
# }

compute_lagged_hfrs_fast <- function(
  hosps,
  deaths,
  lag,
  w = 1,
  real_time = TRUE,
  tvar_lag = NULL,
  dates = NULL
) {
  ddates <- dates %||% names(deaths)
  ddates <- as.Date(ddates)
  deaths <- enframe(deaths, "time_value", "deaths") |>
    mutate(time_value = as.Date(time_value))
  if (!is_null(tvar_lag)) {
    if (length(tvar_lag) != (nhosp <- length(hosps))) {
      cli_abort(
        "`tvar_delay` must have {nhosp} rows. It has {length(tvar_lag)}."
      )
    }
  }
  hosps <- enframe(hosps, "time_value", "hosps") |>
    mutate(time_value = as.Date(time_value), lag = tvar_lag %||% lag)
  if (real_time) {
    # move hosps forward y_t / x_{t-l}, note opposite
    align <- "right"
    hosps <- hosps |>
      mutate(time_value = time_value + lag) |>
      arrange(time_value) |>
      summarise(hosps = mean(hosps), .by = "time_value")
  } else {
    # move deaths back y_{t+l} / x_{t}
    align <- "center"
    deaths <- left_join(
      deaths,
      hosps |> select(time_value, lag),
      by = "time_value"
    ) |>
      mutate(time_value = time_value - lag) |>
      arrange(time_value) |>
      summarise(deaths = mean(deaths), .by = "time_value")
  }

  res <- full_join(hosps, deaths, by = "time_value") |>
    arrange(time_value) |>
    fill(hosps, deaths, .direction = "down") |>
    mutate(
      deaths = roll_mean(deaths, w, align = align, na.rm = TRUE, fill = NA),
      hosps = roll_mean(hosps, w, align = align, na.rm = TRUE, fill = NA),
      HFR = deaths / hosps,
    ) |>
    fill(HFR, .direction = "down") |>
    filter(!is.na(hosps), !is.na(deaths)) |>
    arrange(time_value) |>
    filter(time_value %in% ddates) |>
    select(time_value, HFR) |>
    deframe()
  res
}

tune_lagged_hfrs_retro <- function(
  hosps,
  deaths,
  delay_distr,
  ws,
  lag,
  K = 5,
  tvar_lag = NULL,
  tvar_delay = NULL # needed for CV prediction
) {
  d1 <- ifelse(is_null(tvar_delay), length(delay_distr), ncol(tvar_delay))
  d <- d1 - 1
  death_dates <- as.Date(names(deaths))
  centered <- TRUE

  # --- Phase 1: Fit all HFRs ---
  full_fit <- sapply(ws, function(w) {
    hosps_smoothed <- smooth_counts(hosps, w, centered)
    deaths_smoothed <- smooth_counts(deaths, w, centered)
    lagged_hfrs <- compute_lagged_hfrs_fast(
      hosps_smoothed,
      deaths_smoothed,
      lag,
      dates = death_dates,
      real_time = FALSE,
      tvar_lag = tvar_lag
    )
  })

  # --- Phase 2: Get CV errors ---
  # --- Phase 2a: Preprocess ---
  # Get deaths to hold out
  n_held_out_deaths <- floor(length(death_dates) / K)
  held_out_death_list <- lapply(1:K, function(k) {
    0:(n_held_out_deaths - 1) * K + k
  })
  deaths_CV <- lapply(1:K, function(k) {
    deaths[held_out_death_list[[k]]] <- NA
    deaths_imputed <- round(na.approx(deaths, rule = 2))
    names(deaths_imputed) <- names(deaths)
    return(deaths_imputed)
  })

  # Get corresponding hospitalizations to hold out
  held_out_hosp_list <- lapply(1:K, function(k) {
    new_from = d + k - K * floor((d + k - 1) / K)
    seq(from = new_from, to = length(hosps), by = K)
  })
  hosps_CV <- lapply(1:K, function(k) {
    hosps[held_out_hosp_list[[k]]] <- NA
    hosps_imputed <- round(na.approx(hosps, rule = 2))
    names(hosps_imputed) <- names(hosps)
    return(hosps_imputed)
  })

  # --- Phase 2b: Get K x |ws| matrix of CV errors ---
  cv_errors_list <- lapply(1:K, function(k) {
    # Compute lagged HFRs in each fold
    lagged_estimates <- lapply(ws, function(w) {
      # Smooth the data
      hosps_smoothed <- smooth_counts(hosps_CV[[k]], w, centered)
      deaths_smoothed <- smooth_counts(deaths_CV[[k]], w, centered)
      compute_lagged_hfrs_fast(
        hosps_smoothed,
        deaths_smoothed,
        lag,
        dates = death_dates,
        real_time = FALSE,
        tvar_lag = tvar_lag
      )
    })

    # Compute CV errors on the held-out data ---
    held_out_death_dates <- death_dates[held_out_death_list[[k]]]
    Y_held_out <- deaths[as.Date(names(deaths)) %in% held_out_death_dates]
    # Create a vector of CV errors (by ws)
    cv_errors_k <- sapply(lagged_estimates, function(ests) {
      EY <- compute_expected_deaths(
        hosps,
        ests,
        delay_distr,
        held_out_death_dates,
        tvar_delay = tvar_delay
      )
      mae(EY, Y_held_out)
    })
    cv_errors_k
  })
  cv_errors <- do.call(rbind, cv_errors_list)
  colnames(cv_errors) <- ws

  results <- list(
    full_fit = full_fit,
    cv_errors = cv_errors,
    cvstats = cverr(cv_errors, ws)
  )
  return(results)
}

tune_lagged_hfrs_rt <- function(
  hosps,
  deaths,
  t,
  delay_distr,
  ws,
  lag,
  nstore = 7,
  tvar_lag = NULL,
  tvar_delay = NULL # needed for CV prediction
) {
  # Tune via forward validation
  d1 <- ifelse(is_null(tvar_delay), length(delay_distr), ncol(tvar_delay))
  d <- d1 - 1
  death_dates <- as.Date(names(deaths))
  centered <- FALSE

  # --- Phase 1: Compute HFR estimates ---
  full_fit <- sapply(ws, function(w) {
    # Smooth the data
    hosps_smoothed <- smooth_counts(hosps, w, centered)
    deaths_smoothed <- smooth_counts(deaths, w, centered)

    compute_lagged_hfrs_fast(
      hosps_smoothed,
      deaths_smoothed,
      lag,
      real_time = TRUE,
      tvar_lag = tvar_lag
    )
  })

  # --- Phase 2: Compute FV errors and aggregates ---
  # dims: dates x ws
  n <- nrow(full_fit)
  recent <- (n - 7 + 1):n
  imputed_hfrs <- full_fit[-recent, ]
  slopes <- diff(tail(imputed_hfrs, 2))
  imputed_hfrs <- rbind(
    imputed_hfrs,
    clamp(rep(1, 7) %*% tail(imputed_hfrs, 1) + 1:7 %*% slopes)
  )

  if (is_null(tvar_delay)) {
    Zmat <- make_Z_mat_simple(n, hosps, delay_distr)[recent, ]
  } else {
    delay_mat <- make_tvar_delay_mat(tvar_delay, submat = TRUE)
    Zmat <- make_Z_mat_new(hosps, delay_mat, recent)
  }
  EY <- Zmat[, -(1:d)] %*% imputed_hfrs
  fv_errors_mat <- abs(tail(deaths, 7) - EY)

  list(
    full_fit = tail(full_fit, nstore),
    fv_errors = fv_errors_mat,
    fvstats = cverr(fv_errors_mat, ws)
  )
}


# Convolutional ratio ---------------------------------------------------------

compute_conv_hfr_realtime <- function(hosps, deaths, delay_distr, t) {
  d <- length(delay_distr) - 1
  hosps_in_trailing_window <- hosps[names(hosps) %in% seq(t - d, t, by = "day")]
  if (length(hosps_in_trailing_window) < d + 1) {
    return(NA)
  }
  contributing_hosps <- sum(rev(hosps_in_trailing_window) * delay_distr)

  deaths_at_t <- sum(deaths[(names(deaths) == t)])
  deaths_at_t / contributing_hosps
}

compute_conv_hfrs_realtime <- function(
  hosps,
  deaths,
  delay_distr,
  dates,
  tvar_delay = NULL
) {
  nhosps <- length(hosps)
  ndeaths <- length(deaths)
  d1 <- ifelse(is_null(tvar_delay), length(delay_distr), ncol(tvar_delay))
  d <- d1 - 1
  if (is_null(tvar_delay)) {
    stopifnot(ndeaths > d, nhosps > d)
    delay_mat <- make_delay_mat0(ndeaths, nhosps, delay_distr)
  } else {
    if (nrow(tvar_delay) != nhosps) {
      cli_abort(
        "`tvar_delay` must have {nhosps} rows. It has {nrow(tvar_delay)}."
      )
    }
    delay_mat <- make_tvar_delay_mat(tvar_delay, submat = TRUE)
  }
  Zmat <- make_Z_mat_new(hosps, delay_mat)
  conv_hfrs <- deaths / rowSums(Zmat)
  conv_hfrs[is.nan(conv_hfrs)] <- 0
  names(conv_hfrs) <- dates
  return(conv_hfrs)
}

compute_conv_hfrs_retro <- function(
  hosps,
  deaths,
  delay_distr,
  tvar_delay = NULL
) {
  # https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010554
  ndeaths <- length(deaths)
  nhosps <- length(hosps)
  d1 <- ifelse(is_null(tvar_delay), length(delay_distr), ncol(tvar_delay))
  d <- d1 - 1
  nhfrs <- ndeaths - d
  if (is_null(tvar_delay)) {
    stopifnot(ndeaths > d, nhosps > d)
    Zmat <- make_Z_mat_simple(ndeaths, hosps, delay_distr)
    Dmat <- t(make_delay_mat0(ndeaths, ndeaths, delay_distr))[1:nhfrs, ]
  } else {
    if (nrow(tvar_delay) != nhosps) {
      cli_abort(
        "`tvar_delay` must have {nhosps} rows. It has {nrow(tvar_delay)}."
      )
    }
    delay_mat <- make_tvar_delay_mat(tvar_delay, submat = TRUE)
    Zmat <- make_Z_mat_new(hosps, delay_mat)
    Dmat <- spread_tvar_delay(tvar_delay[1:ndeaths, ])
  }
  denoms_after_t <- 1 / rowSums(Zmat)
  num_after_t <- colScale(Dmat, deaths)

  forward_conv_hfr <- drop(num_after_t %*% denoms_after_t)
  names(forward_conv_hfr) <- names(deaths)[1:nhfrs]
  return(forward_conv_hfr)
}

tune_conv_hfrs_retro <- function(
  hosps,
  deaths,
  delay_distr,
  ws,
  K = 5,
  tvar_delay = NULL
) {
  d1 <- ifelse(is_null(tvar_delay), length(delay_distr), ncol(tvar_delay))
  d <- d1 - 1
  death_dates <- as.Date(names(deaths))
  centered <- TRUE

  # --- Phase 1: Compute HFR estimates on all ws ---
  full_fit <- sapply(ws, function(ww) {
    hosps_smoothed <- smooth_counts(hosps, ww, centered)
    deaths_smoothed <- smooth_counts(deaths, ww, centered)
    tuned_hfrs <- compute_conv_hfrs_retro(
      hosps_smoothed,
      deaths_smoothed,
      delay_distr,
      tvar_delay = tvar_delay
    )
    tuned_hfrs
  })

  # --- Phase 2: Compute CV errors ---
  # --- Phase 2a: Preprocess ---
  # Get deaths to hold out
  n_held_out_deaths <- floor(length(death_dates) / K)
  held_out_death_list <- lapply(1:K, function(k) {
    0:(n_held_out_deaths - 1) * K + k
  })
  deaths_CV <- lapply(1:K, function(k) {
    deaths[held_out_death_list[[k]]] <- NA
    deaths_imputed <- round(na.approx(deaths, rule = 2))
    names(deaths_imputed) <- names(deaths)
    return(deaths_imputed)
  })

  # Get corresponding hospitalizations to hold out
  held_out_hosp_list <- lapply(1:K, function(k) {
    new_from = d + k - K * floor((d + k - 1) / K)
    seq(from = new_from, to = length(hosps), by = K)
  })
  hosps_CV <- lapply(1:K, function(k) {
    hosps[held_out_hosp_list[[k]]] <- NA
    hosps_imputed <- round(na.approx(hosps, rule = 2))
    names(hosps_imputed) <- names(hosps)
    return(hosps_imputed)
  })

  # --- Phase 2b: Compute CV error matrix ---
  cv_errors_list <- lapply(1:K, function(k) {
    # --- Compute HFR estimates on in-fold data ---
    conv_estimates <- lapply(ws, function(w) {
      # Smooth the data
      hosps_smoothed <- smooth_counts(hosps_CV[[k]], w, centered)
      deaths_smoothed <- smooth_counts(deaths_CV[[k]], w, centered)
      compute_conv_hfrs_retro(
        hosps_smoothed,
        deaths_smoothed,
        delay_distr,
        tvar_delay = tvar_delay
      )
    })

    # Compute CV errors on out-of-fold data
    held_out_death_dates <- death_dates[held_out_death_list[[k]]]
    Y_held_out <- deaths[as.Date(names(deaths)) %in% held_out_death_dates]
    cv_errors_k <- sapply(conv_estimates, function(conv_hfrs_cv) {
      # NAs at death dates where HFR dates aren't equipped to predict
      EY <- compute_expected_deaths(
        hosps,
        conv_hfrs_cv,
        delay_distr,
        held_out_death_dates,
        tvar_delay = tvar_delay
      )
      mae(EY, Y_held_out)
    })
    cv_errors_k
  })
  cv_errors <- do.call(rbind, cv_errors_list)
  colnames(cv_errors) <- ws

  results <- list(
    full_fit = full_fit,
    cv_errors = cv_errors,
    cvstats = cverr(cv_errors, ws)
  )
  return(results)
}

tune_conv_hfrs_rt <- function(
  hosps,
  deaths,
  t,
  delay_distr,
  ws,
  nstore = 7,
  tvar_delay = NULL
) {
  # Tune via forward validation
  d1 <- ifelse(is_null(tvar_delay), length(delay_distr), ncol(tvar_delay))
  d <- d1 - 1
  hfr_dates_daily <- seq(t - 7 - d, t, by = "day")
  death_dates <- as.Date(names(deaths))
  centered <- FALSE

  # --- Phase 1: Compute HFR estimates ---
  full_fit <- sapply(ws, function(w) {
    # Smooth the data
    hosps_smoothed <- smooth_counts(hosps, w, centered)
    deaths_smoothed <- smooth_counts(deaths, w, centered)

    conv_hfrs_cv <- compute_conv_hfrs_realtime(
      hosps_smoothed,
      deaths_smoothed,
      delay_distr,
      hfr_dates_daily,
      tvar_delay = tvar_delay
    )
    conv_hfrs_cv <- na.approx(conv_hfrs_cv, na.rm = FALSE, rule = 2)
    names(conv_hfrs_cv) <- hfr_dates_daily
    return(conv_hfrs_cv)
  })

  # --- Phase 2: Compute FV errors (matrix) and aggregate ---
  # Done by linearly extrapolating HFR over the last 7 days
  n <- nrow(full_fit)
  recent <- (n - 7 + 1):n
  imputed_hfrs <- full_fit[-recent, ]
  slopes <- diff(tail(imputed_hfrs, 2))
  imputed_hfrs <- rbind(
    imputed_hfrs,
    clamp(rep(1, 7) %*% tail(imputed_hfrs, 1) + 1:7 %*% slopes)
  )
  if (is_null(tvar_delay)) {
    Zmat <- make_Z_mat_simple(n, hosps, delay_distr)[recent, ]
  } else {
    delay_mat <- make_tvar_delay_mat(tvar_delay, submat = TRUE)
    Zmat <- make_Z_mat_new(hosps, delay_mat, recent)
  }
  EY <- Zmat[, -(1:d)] %*% imputed_hfrs
  fv_errors_mat <- abs(tail(deaths, 7) - EY)

  list(
    full_fit = tail(full_fit, nstore),
    fv_errors = fv_errors_mat,
    fvstats = cverr(fv_errors_mat, ws)
  )
}

# Trend filter helper functions -----------------------------------------------

# make_delay_mat <- function(ndeaths, nhosps, delay_distr) {
#   d <- length(delay_distr) - 1
#   stopifnot(ndeaths + d == nhosps)
#   delay_distr <- rev(delay_distr)
#   M <- Matrix(0, nrow = ndeaths, ncol = nhosps, sparse = TRUE)
#   for (i in 1:ndeaths) {
#     M[i, i:(i + d)] <- delay_distr
#   }
#   M
# }

make_delay_mat0 <- function(ndeaths, nhosps, delay_distr) {
  d <- length(delay_distr)
  zerosr <- rep(0, nhosps - d)
  zerosl <- rep(0, ndeaths - 1)
  # x goes top right to bottom left
  drop0(toeplitz2(c(zerosl, delay_distr, zerosr), ndeaths, nhosps))
}

make_Z_mat_new <- function(
  all_hosps,
  delay_mat,
  idr = seq_len(nrow(delay_mat)),
  idc = seq_along(all_hosps)
) {
  delay_mat <- delay_mat[, idc]
  rs <- 1 / rowSums(delay_mat[idr, ])
  dimScale(delay_mat[idr, ], rs, all_hosps[idc])
}

make_Z_mat_simple <- function(ndeaths, hosps, delay_distr) {
  delay_mat <- make_delay_mat0(ndeaths, length(hosps), delay_distr)
  make_Z_mat_new(hosps, delay_mat)
}

tail_regularizer_mat <- function(n, delay, locs = seq_len(n)) {
  full_delay <- c(delay, rep(0, n - length(delay)))
  partial_delay <- rev(rev(full_delay)[locs])
  partial_delay <- partial_delay / sum(partial_delay)
  pos_vals <- partial_delay[partial_delay > 0]
  partial_delay[partial_delay > 0] <- 1 / sqrt(cumsum(pos_vals))
  drop0(Diagonal(x = rev(partial_delay)[-1]))
}

make_tvar_delay_mat <- function(delay_mat, submat = TRUE) {
  # delay_mat (pi) is an (n x d+1) matrix, P[i,k] = p_{i,k}
  n <- nrow(delay_mat) # nhosps
  d1 <- ncol(delay_mat) # d + 1
  d <- d1 - 1
  R <- t(bandSparse(n, n, k = 0:d, diagonals = delay_mat))
  # R <- Matrix(0, n, n)
  # for (i in 1:n) {
  #   for (j in 1:i) {
  #     if (i - j >= 0 && i - j < d1) {
  #       R[i, j] <- delay_mat[j, i - j + 1]
  #     }
  #   }
  # }
  if (submat) {
    R <- R[d1:n, ]
  }
  R
}

spread_tvar_delay <- function(delay_mat) {
  n <- nrow(delay_mat)
  p <- ncol(delay_mat)
  d <- p - 1
  bandSparse(n, n + d, 0:d, diagonals = delay_mat)[-(1:d), -(1:d)]
  # Dmat <- Matrix(0, n - d, n)
  # for (iter in 1:(n - d)) {
  #   ridx <- iter + d
  #   Dmat[iter, iter:ridx] <- delay_mat[ridx, ]
  # }
  # Dmat
}


compute_lambda_max <- function(
  Z,
  Y,
  order,
  family = c("poissonLI", "gaussian"),
  approx = TRUE
) {
  family <- match.arg(family)
  Dmat <- get_diff_mat(ncol(Z), order)
  qrM <- qr(t(Dmat))
  if (approx) {
    if (family == "poissonLI") {
      ones <- rep(1, nrow(Z))
      return(max(abs(solve(qrM, crossprod(Z, ones)))) * 100)
    } else {
      return(max(abs(solve(qrM, crossprod(Z, Y)))))
    }
  } else {
    stop("Exact lambda_max not implemented yet")
  }
}


# Trend filtering -------------------------------------------------------------

compute_tf_hfrs <- function(
  Z,
  Y,
  lambda,
  order = 2,
  poisson = TRUE,
  var_method = "Estimated",
  deaths_var = NULL,
  gamma = 0,
  real_time = FALSE,
  d = length(delay_distr_T) - 1,
  natural_order = NULL,
  delay_distr_T,
  hosp_locs = seq_len(ncol(Z))
) {
  n_hfrs_until_t <- ncol(Z)
  diff_mat <- get_diff_mat(n_hfrs_until_t, order + 1, hosp_locs)
  death_dates_until_t <- as.Date(names(Y))

  Thetas <- Variable(n_hfrs_until_t)
  if (real_time) {
    # Extra regularization: Towards flat tail. Prepare matrix
    diff_mat_one <- get_diff_mat(n_hfrs_until_t, 1, hosp_locs)
    W <- tail_regularizer_mat(n_hfrs_until_t, delay_distr_T, hosp_locs)
  }

  if (poisson) {
    constraints <- list(Thetas >= 0)
    loss <- mean(
      as_cvxr_expr(Z) %*% Thetas - Y * log(as_cvxr_expr(Z) %*% Thetas)
    )
    tf_regularizer <- lambda * mean(abs(as_cvxr_expr(diff_mat) %*% Thetas))
    loss <- loss + tf_regularizer

    if (real_time) {
      taper_errors <- (W %*% as_cvxr_expr(diff_mat_one) %*% Thetas)**2
      taper_regularizer <- gamma * mean(taper_errors)
      loss <- loss + taper_regularizer
    }
  } else {
    if (var_method == "Estimated") {
      vars <- seven_day_smoothing(Y, centered = TRUE)
    } else if (var_method == "Oracle") {
      vars <- deaths_var[as.Date(names(deaths_var)) %in% death_dates_until_t]
    } else {
      vars <- rep(mean(Y), length(Y))
    }
    if (min(vars) < 1) {
      # Shift upwards so minimum is 1
      vars <- vars + (1 - min(vars))
    }

    z <- Variable(nrow(diff_mat)) # n_hfrs_until_t - order + 1
    constraints <- list(
      z >= 0,
      Thetas >= 0,
      as_cvxr_expr(diff_mat) %*% Thetas >= -z,
      as_cvxr_expr(diff_mat) %*% Thetas <= z
    )
    #   loss <- mean((1/vars)*(Y - Z %*% Thetas)**2) + lambda*mean(z)
    # }
    loss <- mean((1 / vars) * (Y - as_cvxr_expr(Z) %*% Thetas)**2)
    tf_regularizer <- lambda * mean(z)
    loss <- loss + tf_regularizer

    if (real_time) {
      q <- Variable(nrow(W))
      constraints <- append(
        constraints,
        c(
          q >= 0,
          as_cvxr_expr(W) %*% as_cvxr_expr(diff_mat_one) %*% Thetas >= -q,
          as_cvxr_expr(W) %*% as_cvxr_expr(diff_mat_one) %*% Thetas <= q
        )
      )
      taper_regularizer <- gamma * mean(q)
      loss <- loss + taper_regularizer
    }
  }

  # Extra regularization: Natural trend filtering
  if (real_time) {
    # No knots at tail, unless pre-specified (had been 1)
    tail_constraint_order <- ifelse(
      is.null(natural_order),
      order,
      natural_order
    )
    ncd <- ncol(diff_mat)
    last_dates <- Matrix(0, nrow = 1, ncol = ncd, sparse = TRUE)
    if (tail_constraint_order == 0) {
      # Flat at tail
      last_dates[, ncd - 1:0] <- c(1, -1)
    }
    if (tail_constraint_order == 1) {
      # Linear at tail
      last_dates[, ncd - 2:0] <- c(1, -2, 1)
    }
    if (tail_constraint_order == 2) {
      # Quadratic at tail
      last_dates[, ncd - 3:0] <- c(1, -3, 3, -1)
    }
    constraints <- append(constraints, as_cvxr_expr(last_dates) %*% Thetas == 0)
  }

  prob <- Problem(Minimize(loss), constraints = constraints)
  sol <- psolve(prob, "CLARABEL", verbose = FALSE) # Still writes, oh well

  if (is.null(sol)) {
    return(rep(NA, n_hfrs_until_t))
  }
  if (status(prob) != "optimal") {
    print(status(prob))
  }
  if (status(prob) == "solver_error") {
    return(rep(NA, n_hfrs_until_t))
  }
  tf_hfrs <- drop(value(Thetas))
  return(tf_hfrs)
}

# compared to old tune_lam_tf()
# ~~1. this remakes the Z matrix on the fly~~
# 2. keeps the cv error matrix
# 3. returns the fit at all lambdas
# 4. returns summaries of the CV path
cv_lambda_tf <- function(
  hosps,
  deaths,
  delay_distr,
  tvar_delay = NULL,
  lambdas = NULL,
  lambda_max = NULL,
  lambda_min_ratio = 1e-6,
  rescale = TRUE,
  gamma = 0,
  nlambdas = 30,
  K = 5,
  order = 2,
  real_time = FALSE,
  var_method = "Estimated",
  deaths_var = NULL,
  natural_order = order,
  poisson = TRUE
) {
  ndeath <- length(deaths)
  nhosp <- length(hosps)
  d1 <- ifelse(is_null(tvar_delay), length(delay_distr), ncol(tvar_delay))
  d <- d1 - 1

  if (is_null(tvar_delay)) {
    stopifnot(ndeath + d == nhosp)
    Z <- make_Z_mat_simple(ndeath, hosps, delay_distr)
    delay_distr_T <- delay_distr
  } else {
    if (nrow(tvar_delay) != nhosp) {
      cli_abort(
        "`tvar_delay` must have {nhosp} rows. It has {nrow(tvar_delay)}."
      )
    }
    delay_mat <- make_tvar_delay_mat(tvar_delay, submat = TRUE)
    Z <- make_Z_mat_new(hosps, delay_mat)
    delay_distr_T <- drop(tail(tvar_delay, 1))
  }
  Y <- deaths
  if (rescale) {
    maxmax <- max(max(abs(Z)), max(abs(Y)))
    Z <- Z / maxmax
    Y <- Y / maxmax
  }
  if (is.null(lambdas)) {
    if (is.null(lambda_max)) {
      lambda_max <- compute_lambda_max(
        Z,
        Y,
        order,
        family = ifelse(poisson, "poissonLI", "gaussian")
      )
    }
    lambdas <- geom_spacing(lambda_max * lambda_min_ratio, lambda_max, nlambdas)
  }

  full_fit <- sapply(lambdas, \(lambda) {
    compute_tf_hfrs(
      Z,
      Y,
      lambda = lambda,
      order = order,
      poisson = poisson,
      real_time = real_time,
      gamma = gamma,
      d = d,
      var_method = var_method,
      deaths_var = deaths_var,
      natural_order = natural_order,
      delay_distr_T = delay_distr_T
    )
  })

  # never omit the tails, in case of real-time
  middle_fold <- rep(1:K, length.out = ndeath - 2 * (order + 2))
  fold_deaths <- c(rep(0, order + 2), middle_fold, rep(0, order + 2))
  fold_hosps <- c(rep(0, d), fold_deaths)
  cv_errors <- matrix(nrow = K, ncol = length(lambdas))
  for (k in seq_len(K)) {
    train_deaths <- fold_deaths != k
    # train_hosps <- fold_hosps != k
    Ycv <- Y[train_deaths]
    # xcv <- hosp_locs[train_hosps]
    # Zcv <- make_Z_mat(hosps, delay_mat, idr = which(train_deaths), idc = xcv)
    # if (rescale) Zcv <- Zcv / maxmax
    cv_error <- map_dbl(lambdas, \(lambda) {
      tf_hfrs <- tryCatch(
        {
          # Setting real-time=TRUE runs natural TF without tapered smoothing
          compute_tf_hfrs(
            Z[train_deaths, ],
            Ycv,
            lambda = lambda,
            order = order,
            poisson = poisson,
            var_method = var_method,
            deaths_var = deaths_var,
            real_time = FALSE,
            d = d,
            natural_order = natural_order,
            delay_distr_T = delay_distr_T,
            gamma = gamma
          )
        },
        error = function(cond) {
          return(NA)
        }
      )
      if (is.null(tf_hfrs)) {
        return(NA)
      }
      Y_pred <- Z[!train_deaths, ] %*% tf_hfrs
      return(mae(Y_pred, Y[!train_deaths]))
    })
    cv_errors[k, ] <- cv_error
  }
  cvstats <- cverr(cv_errors, lambdas)
  return(list(
    full_fit = full_fit,
    cv_errors = cv_errors,
    K = K,
    lambdas = lambdas,
    cvstats = cvstats
  ))
}

#' takes in a matrix of cv errors, 1 column per lambda, 1 row per fold
#' computes minimizer and 1se lambdas, as well as their indexes
#'
cverr <- function(cvs, lambda) {
  cvm <- colMeans(cvs, na.rm = TRUE)
  cvse <- apply(cvs, 2, sd, na.rm = TRUE) / sqrt(nrow(cvs))
  i0 = which.min(cvm)
  if (length(i0) == 0) {
    i0 <- max(lambda)
  }
  lam.min = lambda[i0]
  # lam.1se = max(lambda[cvm <= cvm[i0] + cvse[i0]], max(lambda), na.rm = TRUE)
  cand <- lambda[cvm <= cvm[i0] + cvse[i0]]
  lam.1se <- if (length(cand) > 0 && any(!is.na(cand))) {
    max(cand, na.rm = TRUE)
  } else {
    max(lambda, na.rm = TRUE)
  }
  i.min = which(lambda == lam.min)
  i.1se = which(lambda == lam.1se)
  return(structure(
    list(
      cvm = cvm,
      cvse = cvse,
      up = cvm + cvse,
      down = pmax(cvm - cvse, 0, na.rm = TRUE),
      lam.min = lam.min,
      lam.1se = lam.1se,
      i.min = i.min,
      i.1se = i.1se,
      lambdas = lambda
    ),
    class = "cvstats"
  ))
}


plot.cvstats <- function(x, ...) {
  rng <- range(x$up, x$down)
  plot(
    x = x$lambdas,
    y = x$cvm,
    pch = 16,
    col = "darkblue",
    log = "x",
    ylim = rng,
    ...
  )
  segments(x$lambdas, x$down, x$lambdas, x$up, col = "darkred")
  abline(v = c(x$lam.min, x$lam.1se), lty = "dashed")
}

fv_gamma_tf <- function(
  hosps,
  deaths,
  lambda,
  gammas,
  delay_distr,
  tvar_delay = NULL,
  window_length = 7,
  order = 2,
  poisson = TRUE,
  var_method = "Estimated",
  deaths_var = NULL,
  natural_order = order,
  rescale = 1
) {
  ndeath <- length(deaths)
  nhosp <- length(hosps)
  d1 <- ifelse(is_null(tvar_delay), length(delay_distr), ncol(tvar_delay))
  d <- d1 - 1

  if (is_null(tvar_delay)) {
    stopifnot(ndeath + d == nhosp)
    Z <- make_Z_mat_simple(ndeath, hosps, delay_distr) / rescale
    delay_distr_T <- delay_distr
  } else {
    if (nrow(tvar_delay) != nhosp) {
      cli_abort(
        "`tvar_delay` must have {nhosp} rows. It has {nrow(tvar_delay)}."
      )
    }
    delay_mat <- make_tvar_delay_mat(tvar_delay, submat = TRUE)
    Z <- make_Z_mat_new(hosps, delay_mat) / rescale
    delay_distr_T <- drop(tail(tvar_delay, 1))
  }
  Y <- deaths / rescale
  n_death_dates <- length(Y)

  errors_all_gamma <- matrix(NA, nrow = window_length, ncol = length(gammas))
  for (days_ago in window_length:1) {
    n_deaths_to_s <- n_death_dates - days_ago
    n_hfrs_to_s <- n_death_dates + d - days_ago
    rows_to_s <- seq_len(n_deaths_to_s)
    cols_to_s <- seq_len(n_hfrs_to_s)
    Zs <- Z[rows_to_s, cols_to_s]
    Ys <- Y[rows_to_s]

    errors_s1 <- map_dbl(gammas, \(gamma) {
      tf_hfrs_fv <- tryCatch(
        {
          compute_tf_hfrs(
            Zs,
            Ys,
            lambda,
            gamma = gamma,
            order = order,
            poisson = poisson,
            real_time = TRUE,
            d = d,
            var_method = var_method,
            deaths_var = deaths_var,
            natural_order = natural_order,
            delay_distr_T = delay_distr
          )
        },
        error = function(cond) {
          return(NA)
        }
      )
      if (is.null(tf_hfrs_fv)) {
        return(NA)
      }
      recent_hfrs <- tf_hfrs_fv[(n_hfrs_to_s - d + 1):n_hfrs_to_s]
      imputed_hfr <- recent_hfrs[d] + (recent_hfrs[d] - recent_hfrs[d - 1])
      recent_hfrs <- append(recent_hfrs, imputed_hfr)

      # Estimate deaths at s+1 with these HFR estimates
      hosp_conv_s1 <- Z[
        (n_deaths_to_s + 1),
        (n_hfrs_to_s - d + 1):(n_hfrs_to_s + 1)
      ]
      Y_reconvolved <- sum(hosp_conv_s1 * recent_hfrs)

      # Take absolute error
      Y_s1 <- Y[n_deaths_to_s + 1]
      return(abs(Y_reconvolved - Y_s1))
    })
    errors_all_gamma[window_length - days_ago + 1, ] <- errors_s1
  }
  errors_all_gamma
}

tune_trendfilter_rt <- function(
  hosps,
  deaths,
  gammas,
  delay_distr,
  tvar_delay = NULL,
  lambdas = NULL,
  order = 2,
  poisson = TRUE,
  K = 5,
  which_lambdas = c("min", "1se"),
  nlambdas = 20,
  lambda_max = NULL,
  lambda_min_ratio = 1e-4,
  rescale = TRUE,
  var_method = "Estimated",
  deaths_var = NULL,
  natural_order = order,
  fv_errors = NULL,
  window_length = 7,
  nstore = 7
) {
  which_lambdas <- match.arg(which_lambdas, several.ok = TRUE)
  if (length(lambdas) == 1) {
    which_lambdas <- "fixed"
  } else {
    # Tune lambda
    tune_lambda_gamma0 <- cv_lambda_tf(
      hosps,
      deaths,
      delay_distr,
      tvar_delay = tvar_delay,
      lambdas = lambdas,
      lambda_max = lambda_max,
      lambda_min_ratio = lambda_min_ratio,
      rescale = rescale,
      var_method = var_method,
      deaths_var = deaths_var,
      nlambdas = nlambdas,
      K = K,
      order = order,
      natural_order = natural_order,
      poisson = poisson,
      real_time = FALSE
    )
    lambdas <- lambdas %||% tune_lambda_gamma0$lambdas
  }
  if (is_null(tvar_delay)) {
    Z <- make_Z_mat_simple(length(deaths), hosps, delay_distr)
    delay_distr_T <- delay_distr
  } else {
    delay_mat <- make_tvar_delay_mat(tvar_delay, submat = TRUE)
    Z <- make_Z_mat_new(hosps, delay_mat)
    delay_distr_T <- drop(tail(tvar_delay, 1))
  }
  Y <- deaths
  if (rescale) {
    maxmax <- max(max(abs(Z)), max(abs(Y)))
    Z <- Z / maxmax
    Y <- Y / maxmax
  }
  fv_partial <- function(chosen_lambda) {
    fv_gamma_tf(
      hosps,
      deaths,
      chosen_lambda,
      gammas,
      delay_distr,
      tvar_delay = tvar_delay,
      window_length = window_length,
      order = order,
      poisson = poisson,
      var_method = var_method,
      deaths_var = deaths_var,
      natural_order = natural_order,
      rescale = ifelse(rescale, maxmax, 1)
    )
  }
  which_lambda_vals <- map_dbl(which_lambdas, \(x) {
    switch(
      x,
      fixed = lambdas,
      min = tune_lambda_gamma0$cvstats$lam.min,
      `1se` = tune_lambda_gamma0$cvstats$lam.1se
    )
  }) |>
    set_names(which_lambdas)

  fv_errors_window <- map(which_lambda_vals, \(lambda) fv_partial(lambda))
  fv_stats_window <- map(fv_errors_window, \(errors) cverr(errors, gammas))

  full_fit <- map(unique(which_lambda_vals), \(lambda) {
    sapply(gammas, \(gamma) {
      out <- tryCatch(
        {
          compute_tf_hfrs(
            Z,
            Y,
            lambda = lambda,
            gamma = gamma,
            order = order,
            poisson = poisson,
            real_time = TRUE,
            d = d,
            var_method = var_method,
            deaths_var = deaths_var,
            natural_order = natural_order,
            delay_distr_T = delay_distr_T
          )
        },
        error = function(cond) {
          return(rep(NA, ncol(Z)))
        }
      )
      if (is.null(out)) {
        out <- rep(NA, ncol(Z))
      }
      tail(out, nstore)
    })
  })
  if (length(which_lambdas) != length(full_fit)) {
    full_fit <- list("min" = full_fit[[1]], "1se" = full_fit[[1]])
  } else {
    names(full_fit) <- which_lambdas
  }
  output <- list(
    full_fit = full_fit,
    fv_errors = fv_errors_window,
    fv_stats_window = fv_stats_window
  )
  if (length(which_lambdas) == 1) {
    output <- map(output, 1)
  }
  output
}
