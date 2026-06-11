`%nin%` <- function(x, table) !(x %in% table)

# Real time -------------------------------------------------------------------
find_window_fv <- function(
  res,
  tune_par,
  cvopt = c("none", "min", "1se"),
  window_weeks = 4L
) {
  out <- list()
  cvopt <- match.arg(cvopt)
  errs <- map(res$slide_value, "fv_errors")
  if (cvopt != "none") {
    errs <- map(errs, cvopt)
  }
  for (ii in seq_along(errs)) {
    err_mat <- do.call(rbind, errs[max(1, ii - window_weeks + 1):ii])
    out[[ii]] <- cverr(err_mat, tune_par)
  }
  tibble(
    idx_min = map_dbl(out, "i.min"),
    err_min = map_dbl(out, \(x) x$cvm[x$i.min]),
    idx_1se = map_dbl(out, "i.1se"),
    err_1se = map_dbl(out, \(x) x$cvm[x$i.1se])
  )
}


replace_inf <- function(x, replace = NA) {
  if (any(is.infinite(x))) {
    bad <- is.infinite(x)
    x <- vctrs::vec_assign(x, bad, replace)
  }
  x
}

examine_rt_tuning <- function(
  res,
  job,
  cvopt = c("none", "min", "1se"),
  estimator = c("tf", "lagged", "conv")
) {
  cvopt <- match.arg(cvopt)
  estimator <- match.arg(estimator)
  if (estimator == "tf") {
    tp <- geom_spacing(1, 1e5, 20L)
  } else {
    res <- res[[estimator]]
    cvopt <- "none"
    tp <- c(1, 7, 15, 21, 29)
  }
  tune <- find_window_fv(res, tp, cvopt)
  tune
}

grab_rt_maes <- function(
  res,
  job,
  cvopt = c("none", "min", "1se"),
  fvopt = c("min", "1se"),
  estimator = c("tf", "lagged", "conv"),
  trim_left_tail = 3L
) {
  region <- job$prob.pars$region
  hfrs <- filter(job$problem$data$gt_hfrs, geo_value == region) |>
    select(-geo_value)
  cvopt <- match.arg(cvopt)
  fvopt <- match.arg(fvopt)
  estimator <- match.arg(estimator)
  if (estimator == "tf") {
    tp <- geom_spacing(1, 1e5, 20L)
  } else {
    res <- res[[estimator]]
    cvopt <- "none"
    tp <- c(1, 7, 15, 21, 29)
  }
  idx <- paste("idx", fvopt, sep = "_")
  tune <- find_window_fv(res, tp, cvopt)
  full_fits <- map(res$slide_value, "full_fit")
  if (estimator == "tf") {
    full_fits <- map(full_fits, cvopt)
    which_idx <- tune[[idx]]
  } else {
    which_idx <- pmax(2, tune[[idx]]) # disallow w = 1
  }
  fits <- tibble(
    time_value = res$version,
    est = map2_dbl(full_fits, which_idx, \(x, y) tail(x[, y], 1))
  ) |>
    left_join(hfrs, by = "time_value") |>
    fill(hfr, .direction = "downup")
  if (trim_left_tail > 0) {
    fits <- slice(fits, -(1:trim_left_tail))
  }

  tibble(
    mae = mean(abs(replace_inf(fits$est) - fits$hfr), na.rm = TRUE),
    geo_value = region,
    order = job$algo.pars$order,
    repl = job$repl
  )
}

benchcurve_rt <- function(
  res,
  job,
  fvopt = c("1se", "min"),
  estimator = c("conv", "lagged")
) {
  fvopt <- match.arg(fvopt)
  estimator <- match.arg(estimator)
  res <- res[[estimator]]
  tp <- c(1, 7, 15, 21, 29)
  idx <- paste("idx", fvopt, sep = "_")
  tune <- find_window_fv(res, tp)
  which_idx <- pmax(2, tune[[idx]])
  full_fits <- map(res$slide_value, "full_fit")
  tibble(
    time_value = res$version,
    est = map2_dbl(full_fits, which_idx, \(x, y) tail(x[, y], 1))
  )
}

tune_tf_order <- function(
  res,
  cvopt = c("min", "1se"),
  fvopt = c("min", "1se")
) {
  cvopt <- match.arg(cvopt)
  fvopt <- match.arg(fvopt)
  tp <- geom_spacing(1, 1e5, 20L)
  idx <- paste("idx", fvopt, sep = "_")
  tune <- find_window_fv(res, tp, cvopt)
  full_fits <- map(res$slide_value, \(x) x$full_fit[[cvopt]])
  tibble(
    time_value = res$version,
    est = map2_dbl(full_fits, tune[[idx]], \(x, y) tail(x[, y], 1)),
    cverr = tune$err_min
  )
}

min_order <- function(order, est, cverr) {
  idx <- which(!is.na(est))
  out <- ifelse(any(idx), order[idx][which.min(cverr[idx])], NA)
  out
}


# Retrospective ---------------------------------------------------------------
# In this case, estimates are not always on the same dates, so we must
# first collect curves, then intersect targets.
tfcurve_retro <- function(
  res,
  job,
  which_lambda = c("min", "1se", "oracle")
) {
  region <- job$prob.pars$region
  time_value <- job$problem$data$raw_data |>
    filter(geo_value == region) |>
    pull(time_value)
  which_lambda <- match.arg(which_lambda)
  if (which_lambda %in% c("min", "1se")) {
    idx <- paste("i", which_lambda, sep = ".")
    idx <- res$cvstats[[idx]]
  } else {
    fits <- res$full_fit |> as_tibble(.name_repair = "unique_quiet")
    fits$time_value = time_value
    hfrs <- filter(job$problem$data$gt_hfrs, geo_value == region)
    local_dat <- left_join(fits, hfrs, by = "time_value")
    maes <- colMeans(
      abs(as.matrix(fits |> select(-time_value)) - local_dat$hfr),
      na.rm = TRUE
    )
    idx <- which.min(maes)
  }
  tibble(
    time_value,
    est = res$full_fit[, idx]
  )
}

benchcurve_retro <- function(
  res,
  job,
  which_w = c("1se", "min", "oracle"),
  estimator = c("conv", "lagged")
) {
  estimator <- match.arg(estimator)
  which_w <- match.arg(which_w)
  res <- res[[estimator]]
  if (which_w %in% c("min", "1se")) {
    idx <- paste("i", which_w, sep = ".")
    idx <- pmax(2, res$cvstats[[idx]]) # disallow w = 1
    fits <- res$full_fit[, idx] |>
      enframe("time_value", "est") |>
      mutate(time_value = as.Date(time_value))
    return(fits)
  }
  # oracle
  fits <- res$full_fit |>
    as_tibble(.name_repair = "unique_quiet", rownames = "time_value") |>
    mutate(time_value = as.Date(time_value))
  region <- job$prob.pars$region
  hfrs <- filter(job$problem$data$gt_hfrs, geo_value == region)
  local_dat <- left_join(fits, hfrs, by = "time_value")
  maes <- colMeans(
    abs(as.matrix(fits |> select(-time_value)) - local_dat$hfr),
    na.rm = TRUE
  )
  idx <- which.min(maes)
  tibble(
    time_value = fits$time_value,
    est = res$full_fit[, idx]
  )
}

examine_retro_tuning <- function(
  res,
  job,
  estimator = c("tf", "conv", "lagged")
) {
  if (estimator %in% c("conv", "lagged")) {
    res <- res[[estimator]]
  }
  out <- res$cvstats
  tibble(
    idx_min = out$i.min,
    err_min = out$cvm[out$i.min],
    idx_1se = out$i.1se,
    err_1se = out$cvm[out$i.1se]
  )
}
