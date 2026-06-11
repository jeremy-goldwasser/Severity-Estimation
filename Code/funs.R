compute_max_acf_lag <- function(hosps, deaths, max_lag = 40L) {
  cc <- ccf(hosps, deaths, lag.max = max_lag, plot = FALSE, na.action = na.omit)
  max_correlation_lag <- which.max(cc$acf[1:max_lag])
  oracle_lag <- abs(cc$lag[max_correlation_lag])
  oracle_lag
}

compute_1se_acf_lag <- function(hosps, deaths, max_lag = 40L) {
  cc <- ccf(hosps, deaths, lag.max = max_lag, plot = FALSE, na.action = na.omit)
  # ac1 <- acf(hosps, lag.max = max_lag, plot = FALSE, na.action = na.omit)$acf
  # ac2 <- acf(deaths, lag.max = max_lag, plot = FALSE, na.action = na.omit)$acf
  # se <- sqrt(sum(ac1 * ac2) / cc$n.used)
  se <- sqrt(1 / cc$n.used)
  i0 <- which.max(cc$acf[1:max_lag])
  cand <- which(cc$acf[1:max_lag] >= cc$acf[i0] - se)
  if (length(cand) == 0) {
    return(0)
  }
  oracle_lag <- abs(cc$lag[max(cand)])
  oracle_lag
}

rmultinom1 <- function(size, prob) {
  size <- max(0, size, na.rm = TRUE)
  drop(rmultinom(1, size, prob))
}

fix_weekly_and_negatives <- function(x, zero_cut = 1, back_negatives = 28) {
  n <- length(x)
  xout <- x
  npz <- 0L
  for (ii in seq_along(x)) {
    if (x[ii] == 0) {
      npz <- npz + 1L
    } else if (x[ii] > 0 && npz > zero_cut) {
      ll <- npz + 1L
      back_vec <- rmultinom1(x[ii], rep(1 / ll, ll))
      xout[ii] <- 0
      xout[ii:(ii - npz)] <- xout[ii:(ii - npz)] + back_vec
      npz <- 0L
    } else if (x[ii] < 0) {
      beg <- min(which.max(x > 0), ii)
      probs <- pmax(xout[ii:beg], 1)
      if (length(probs) > 1) {
        probs[1] <- mean(probs[-1])
      }
      if (all(probs <= 0)) {
        x[ii] <- 0
        npz <- npz + 1L
        next
      }
      back_vec <- rmultinom1(abs(x[ii]), rev(probs))
      xout[beg:ii] <- pmax(xout[beg:ii], 0) + back_vec
    }
  }
  xout
}

back_distribute_outliers <- function(x, otlr_idx, resid, window_length = 60L) {
  if (sum(otlr_idx) == 0) {
    return(x)
  }
  xout <- x
  for (ii in which(otlr_idx)) {
    beg <- max(1, ii - window_length)
    idx <- beg:ii
    probs <- pmax(x[idx], 1)
    back_vec <- rmultinom1(resid[ii], probs)
    xout[idx] <- xout[idx] + back_vec
  }
  xout
}

fix_right_tail <- function(x, window_size = 7L) {
  if (max(x, na.rm = TRUE) <= 0) {
    return(rep(0, length(x)))
  }
  nnas <- 0L
  x <- rev(x)
  for (ii in seq_along(x)) {
    if (is.na(x[ii]) || x[ii] == 0) {
      nnas <- nnas + 1L
    } else {
      last_avail <- round(mean(x[ii:(ii + window_size - 1L)], na.rm = TRUE))
      break
    }
  }
  if (nnas > 0L && last_avail > 0) {
    x[1:nnas] <- rmultinom1(nnas * last_avail, rep(1, nnas))
  }
  rev(x)
}

enlist <- function(...) {
  rlang::dots_list(
    ...,
    .homonyms = "error",
    .named = TRUE,
    .check_assign = TRUE
  )
}

process_tibble_to_named_vctrs <- function(tib, d = 60L) {
  hh <- with(tib, set_names(hosps, time_value))
  dd <- with(tib, set_names(deaths, time_value)) |> magrittr::extract(-c(1:d))
  hfr <- tib[["hfr"]]
  if (!is.null(hfr)) {
    names(hfr) <- tib$time_value
  }
  death_dates <- as.Date(names(dd))
  enlist(hh, dd, death_dates, hfr)
}

find_window_fv_1se <- function(
  ll,
  tune_par,
  window_weeks = 4L,
  nm = "fv_errors",
  transp = FALSE
) {
  out <- list()
  ll <- map(ll, \(x) drop(pluck(x, nm)))
  if (transp) {
    ll <- map(ll, t)
  }
  for (ii in seq_along(ll)) {
    out[[ii]] <- cverr(
      do.call(rbind, ll[max(1, ii - window_weeks + 1):ii]),
      tune_par
    )
  }
  tibble(
    idx_1se = map_dbl(out, "i.1se"),
    err_1se = map_dbl(out, \(x) x$cvm[x$i.1se])
  )
}

fix_null_fits <- function(fits, p = 20) {
  if (!is.list(fits)) {
    return(fits)
  }
  nullfits <- map_lgl(fits, is.null)
  if (all(nullfits)) {
    stop("All fits are null.")
  }
  n <- length(fits[[which.min(nullfits)]])
  fits[which(nullfits)] <- map(1:sum(nullfits), ~ rep(NA, n))
  do.call(cbind, fits)
}

detect_outlr_rollmean <- function(
  x = seq_along(y),
  y,
  n = 42,
  detection_multiplier = 4,
  min_radius = "auto",
  log_transform = FALSE
) {
  if (min_radius == "auto") {
    min_radius <- ceiling(mean(pmax(y, 0), na.rm = TRUE))
    min_radius <- max(min_radius, 10)
  }
  min_lower <- ifelse(log_transform, -Inf, 0)
  z <- as_epi_df(tibble(geo_value = 0, time_value = x, y = y, neg = y < 0)) |>
    mutate(y = pmax(0, y))
  if (log_transform) {
    offset <- as.integer(any(z$y == 0))
    min_radius <- log(min_radius)
    z <- mutate(z, y = log(y + offset))
  }
  z <- z |>
    epi_slide(
      fitted = mean(y, na.rm = TRUE),
      .window_size = n,
      .align = "center"
    ) |>
    dplyr::mutate(resid = y - fitted) |>
    epi_slide(
      rsd = sd(resid, na.rm = TRUE),
      .window_size = n,
      .align = "center"
    ) |>
    mutate(
      upper = fitted + pmax(min_radius, detection_multiplier * rsd),
      otlr = y > upper,
      repl = round(fitted * otlr),
      resid = round(resid * otlr)
    )
  z |> select(otlr, repl, resid)
}
