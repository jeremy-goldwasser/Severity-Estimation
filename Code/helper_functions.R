library(lubridate)
library(extraDistr)
library(VGAM)
library(epiprocess)

################ SIMULATING DEATHS ################
sim_pb_deaths <- function(
  death_dates,
  hosps,
  hfrs,
  delay_distr,
  tvar_delay = NULL
) {
  d1 <- ifelse(is_null(tvar_delay), length(delay_distr), ncol(tvar_delay))
  d <- d1 - 1
  nhosps <- length(hosps)
  ndeaths <- length(death_dates)
  stopifnot(nhosps == length(hfrs))
  if (is_null(tvar_delay)) {
    Zmat <- make_Z_mat_simple(nhosps - d, hfrs, delay_distr)
  } else {
    if (nrow(tvar_delay) != nhosps) {
      cli_abort(
        "`tvar_delay` must have {nhosps} rows. It has {nrow(tvar_delay)}."
      )
    }
    delay_mat <- make_tvar_delay_mat(tvar_delay, submat = TRUE)
    Zmat <- make_Z_mat_new(hfrs, delay_mat)
  }
  ptn <- Zmat != 0
  ptn_t <- t(ptn)
  hosp_mat <- colScale(ptn, unname(hosps))

  sim_deaths <- rbinom(nnzero(ptn), t(hosp_mat)@x, t(Zmat)@x) #rowwise extraction
  sim_deaths <- drop(rowSums(t(sparseMatrix(
    i = ptn_t@i,
    p = ptn_t@p,
    x = sim_deaths,
    index1 = FALSE
  ))))

  names(sim_deaths) <- death_dates
  sim_deaths
}

# sim_pb_deaths_old <- function(death_dates, hosps, hfrs, delay_distr) {
#   d <- length(delay_distr) - 1
#   hosp_dates <- names(hosps)
#   hfr_dates <- names(hfrs)
#   sim_deaths <- sapply(death_dates, function(t) {
#     trailing_t <- seq(t - d, t, by = "day")

#     # Match hospitalization and HFR values for the trailing days
#     hosps_s <- hosps[match(trailing_t, hosp_dates)]
#     p_die_t_given_s <- hfrs[match(trailing_t, hfr_dates)] * rev(delay_distr)

#     # Use vectorized rbinom to generate deaths in one call
#     sum(rbinom(length(trailing_t), hosps_s, p_die_t_given_s), na.rm = TRUE)
#   })

#   names(sim_deaths) <- death_dates
#   return(sim_deaths)
# }

sim_normal_deaths <- function(
  death_dates,
  hosps,
  hfrs,
  delay_shape,
  dispersion,
  N = 1
) {
  # Expected number of deaths at t
  hfr_dates <- as.Date(names(hfrs))
  hosp_dates <- as.Date(names(hosps))
  d <- length(delay_shape) - 1
  mean_and_var <- sapply(death_dates, function(t) {
    trailing_t <- seq(t - d, t, by = "day")
    hosps_trailing <- hosps[hosp_dates %in% trailing_t]
    if (length(hfrs) > 1) {
      first_hfrs <- hfrs[hfr_dates %in% trailing_t]
    } else {
      first_hfrs <- rep(hfrs, length(trailing_t))
    }
    delay_distr <- first_hfrs * rev(delay_shape)
    mu_y <- sum(hosps_trailing * delay_distr)
    var_y <- sum(hosps_trailing * delay_distr * (1 - delay_distr))
    c(mu_y, var_y)
  })
  deaths_noiseless <- mean_and_var[1, ]
  vars <- mean_and_var[2, ]
  deaths_list <- lapply(1:N, function(i) {
    e <- rnorm(n = length(death_dates), sd = sqrt(dispersion * vars))
    deaths <- deaths_noiseless + e
    names(deaths) <- death_dates
    deaths[deaths < 0] <- 0
    return(round(deaths))
  })
  if (N == 1) {
    return(deaths_list[[1]])
  } else {
    return(deaths_list)
  }
}

phi2rho <- function(phi) {
  1 / (1 + phi)
}

sim_bb_deaths <- function(
  death_dates,
  hosps,
  hfrs,
  delay_shape,
  dispersion,
  tvar_delay = NULL
) {
  d1 <- ifelse(is_null(tvar_delay), length(delay_shape), ncol(tvar_delay))
  d <- d1 - 1
  nhosps <- length(hosps)
  stopifnot(nhosps == length(hfrs), length(death_dates) == nhosps - d)
  if (is_null(tvar_delay)) {
    Zmat <- make_Z_mat_simple(nhosps - d, hfrs, delay_shape)
  } else {
    if (nrow(tvar_delay) != nhosps) {
      cli_abort(
        "`tvar_delay` must have {nhosps} rows. It has {nrow(tvar_delay)}."
      )
    }
    delay_mat <- make_tvar_delay_mat(tvar_delay, submat = TRUE)
    Zmat <- make_Z_mat_new(hfrs, delay_mat)
  }
  mu_y <- drop(Zmat %*% hosps)
  Zmatm1 <- Zmat
  Zmatm1@x <- 1 - Zmatm1@x
  s2 <- drop((Zmat * Zmatm1) %*% hosps)
  ptn <- Zmat != 0
  n <- drop(ptn %*% hosps)
  mu_h <- mu_y / n
  newvar <- dispersion * s2
  rho <- (newvar / (n * mu_h * (1 - mu_h)) - 1) / (n - 1)
  bb_deaths <- rbetabinom(length(death_dates), n, mu_h, rho = rho)
  names(bb_deaths) <- death_dates
  bb_deaths
}

# sim_bb_deaths_old <- function(
#   death_dates,
#   hosps,
#   hfrs,
#   delay_shape,
#   dispersion,
#   tvar_delay = NULL
# ) {
#   hosp_dates <- as.Date(names(hosps))
#   hfr_dates <- as.Date(names(hfrs))
#   d <- length(delay_shape) - 1
#   bb_deaths <- sapply(death_dates, function(t) {
#     trailing_t <- seq(t - d, t, by = "day")
#     hosps_trailing <- hosps[hosp_dates %in% trailing_t]
#     # mu_y: Expected number of deaths at t
#     first_hfrs <- hfrs[hfr_dates %in% trailing_t]
#     delay_distr <- first_hfrs * rev(delay_shape)
#     mu_y <- sum(hosps_trailing * delay_distr)

#     s2 <- sum(hosps_trailing * delay_distr * (1 - delay_distr))
#     n <- sum(hosps_trailing) # Total number of people who may die at t
#     mu_h <- mu_y / n
#     newvar <- dispersion * s2
#     rho <- (newvar / (n * mu_h * (1 - mu_h)) - 1) / (n - 1)
#     beta_binom_conc <- rbetabinom(1, n, mu_h, rho = rho)
#   })

#   names(bb_deaths) <- death_dates
#   return(bb_deaths)
# }

###### Misc functions for preprocessing, estimation, and evaluation ######
mae <- function(vec1, vec2) {
  mean(abs(as.numeric(vec1) - as.numeric(vec2)), na.rm = T)
}

geom_spacing <- function(Lmin, Lmax, num) {
  r <- (Lmax / Lmin)^(1 / (num - 1))
  Lmin * r^(1:num - 1)
}

Seq = function(a, b) {
  if (a <= b) {
    return(a:b)
  } else {
    return(integer(0))
  }
}


get_diff_mat <- function(p, k, x = seq_len(p)) {
  dspline::d_mat(k, x)
}

seven_day_smoothing <- function(x, centered = TRUE) {
  smooth_counts(x, 7, centered = centered)
}


# seven_day_smoothing_old <- function(data, centered = TRUE) {
#   n <- length(data)
#   if (centered == FALSE) {
#     smoothed <- sapply(1:n, function(i) {
#       # Will be flat for the first 7 days
#       first <- ifelse(i - 6 < 1, 1, i - 6)
#       last <- ifelse(i - 6 < 1, 7, i)
#       mean(data[first:last])
#     })
#   } else {
#     smoothed <- sapply(1:n, function(i) {
#       # Later, should account for boundary effects
#       first <- ifelse(i - 3 < 1, 1, i - 3)
#       last <- ifelse(i + 3 > n, n, i + 3)
#       mean(data[first:last])
#     })
#   }

#   if (!is.null(names(data))) {
#     names(smoothed) <- names(data)
#   }
#   return(smoothed)
# }

# compute_expected_deaths_old <- function(hosps, hfrs, delay_list, death_dates) {
#   # d <- length(delay_list[[1]]) - 1
#   d <- if (is.list(delay_list)) {
#     length(delay_list[[1]]) - 1
#   } else {
#     length(delay_list) - 1
#   }
#   expected_deaths <- sapply(death_dates, function(t) {
#     EY_t <- sum(sapply(0:d, function(k) {
#       tmk <- as.character(t - k)
#       # delay_distr_tmk <- delay_list[[tmk]]
#       delay_distr_tmk <- if (is.list(delay_list)) {
#         delay_list[[tmk]]
#       } else {
#         delay_list
#       }
#       delay_tmk_to_t <- delay_distr_tmk[k + 1]
#       p_tmk_to_t <- delay_tmk_to_t * hfrs[tmk]
#       X_tmk <- hosps[tmk]
#       EY_t_from_tmk <- X_tmk * p_tmk_to_t
#       # This isn't a list, but sapply returning one. Weird.
#       return(EY_t_from_tmk)
#     }))
#   })
#   names(expected_deaths) <- death_dates
#   return(expected_deaths)
# }

compute_expected_deaths <- function(
  hosps,
  hfrs,
  delay_distr,
  death_dates,
  tvar_delay = NULL
) {
  d1 <- ifelse(is_null(tvar_delay), length(delay_distr), ncol(tvar_delay))
  d <- d1 - 1
  nhosps <- length(hosps)
  ndeaths <- nhosps - d
  if (is_null(tvar_delay)) {
    Zmat <- make_Z_mat_simple(ndeaths, hosps, delay_distr)
  } else {
    if (nrow(tvar_delay) != nhosps) {
      cli_abort(
        "`tvar_delay` must have {nhosps} rows. It has {nrow(tvar_delay)}."
      )
    }
    delay_mat <- make_tvar_delay_mat(tvar_delay, submat = TRUE)
    Zmat <- make_Z_mat_new(hosps, delay_mat)
  }
  ddates <- names(hosps)[-(1:d)]
  keep_rows <- ddates %in% death_dates
  keep_cols <- names(hosps) %in% names(hfrs)
  expected_deaths <- drop(Zmat[keep_rows, keep_cols] %*% hfrs)
  names(expected_deaths) <- death_dates
  expected_deaths
}

compute_oracle_vars <- function(hosps, hfrs, delay_list, death_dates) {
  d <- if (is.list(delay_list)) {
    length(delay_list[[1]]) - 1
  } else {
    length(delay_list) - 1
  }
  vars <- sapply(death_dates, function(t) {
    Var_t <- sum(sapply(0:d, function(k) {
      tmk <- as.character(t - k)
      delay_distr_tmk <- if (is.list(delay_list)) {
        delay_list[[tmk]]
      } else {
        delay_list
      }
      delay_tmk_to_t <- delay_distr_tmk[k + 1]
      p_tmk_to_t <- delay_tmk_to_t * hfrs[tmk]
      X_tmk <- hosps[tmk]
      var_t_from_tmk <- X_tmk * p_tmk_to_t * (1 - p_tmk_to_t)
      return(var_t_from_tmk)
    }))
  })
  names(vars) <- death_dates
  return(vars)
}


######### NEWER ADDITIONS
intersect_dates <- function(dates1, dates2) {
  return(as.Date(intersect(dates1, dates2), origin = "1970-01-01"))
}

# compute_optimal_lag <- function(hosps, deaths, verbose = FALSE) {
#   hosp_dates <- names(hosps)
#   death_dates <- names(deaths)
#   intersect_dates <- as.Date(intersect(hosp_dates, death_dates))
#   deaths_both <- deaths[death_dates %in% intersect_dates]
#   hhs_hosps_both <- hosps[hosp_dates %in% intersect_dates]
#   cc <- ccf(hhs_hosps_both, deaths_both, lag.max = 40, plot = FALSE)
#   max_correlation_lag <- which.max(cc$acf[1:41]) # don't allow leads or negative correlation
#   oracle_lag <- abs(cc$lag[max_correlation_lag])
#   if (verbose) {
#     print(paste0("Correlation-maximizing lag is ", oracle_lag, " days."))
#   }
#   return(oracle_lag)
# }

get_gamma_params <- function(Mean, Sd) {
  Var <- Sd**2
  shape <- (Mean**2) / Var
  rate <- Mean / Var
  return(c(shape, rate))
}

# The below corrects for double interval censoring
# See: https://nfidd.github.io/sismid/sessions/using-delay-distributions-to-model-the-data-generating-process-of-an-epidemic.html#discretising-a-delay-distribution
# The difference is subtracting pgamma(x-1) rather than pgamma(x) as done
# in ddgamma(x).
# Compare:
# plot(make_delay_distr(8, 7, 20))
# points(make_delay_distr_old(8, 7, 20), col = 2)
# make_delay_distr_new <- function(Mean, Sd, d, x = 0:d) {
#   params <- get_gamma_params(Mean = Mean, Sd = Sd)
#   shape <- params[1]
#   rate <- params[2]
#   # DelayShape <- ddgamma(0:d, shape, rate)
#   DelayShape <- pgamma(x + 1, shape, rate) - pgamma(x - 1, shape, rate)
#   DelayShape <- DelayShape / sum(DelayShape)
#   return(DelayShape)
# }

make_delay_distr <- function(Mean, Sd, d) {
  params <- get_gamma_params(Mean = Mean, Sd = Sd)
  shape <- params[1]
  rate <- params[2]
  DelayShape <- ddgamma(0:d, shape, rate)
  DelayShape <- DelayShape / sum(DelayShape)
  return(DelayShape)
}


########### Was in preprocess_for_sim.R ############

rm_outliers <- function(y, detection_multiplier = 3, window_size = 21) {
  # detection_multiplier had been 10
  n_death_dates <- length(y)
  if (window_size %% 2 == 0) {
    n_back <- (window_size - 1) / 2
    n_forward <- (window_size - 1) / 2
  } else {
    n_back <- (window_size / 2) - 1
    n_forward <- window_size / 2
  }
  rolling_median <- sapply(1:n_death_dates, function(i) {
    median(y[max(1, i - n_back):min(n_death_dates, i + n_forward)])
  })
  # plot(death_dates, rolling_median, type="l")
  resids <- y - rolling_median
  q3 <- quantile(resids, .75)
  q1 <- quantile(resids, .25)
  iqr <- q3 - q1
  upper_limit <- rolling_median + detection_multiplier * iqr
  lower_limit <- rolling_median - detection_multiplier * iqr

  y_imp <- y
  y_imp[y > upper_limit] <- upper_limit[y > upper_limit]
  y_imp[y < lower_limit] <- lower_limit[y < lower_limit]
  names(y_imp) <- names(y)
  return(y_imp)
}


rm_outliers_rolling <- function(deaths, detection_multiplier, n = 21) {
  roll_iqr = function(
    z,
    detection_multiplier,
    n = 21,
    min_radius = 0,
    replacement_multiplier = 0,
    min_lower = -Inf
  ) {
    if (typeof(z$y) == "integer") {
      as_type = as.integer
    } else {
      as_type = as.numeric
    }

    epiprocess::epi_slide(
      z,
      roll_iqr = stats::IQR(resid),
      before = floor((n - 1) / 2),
      after = ceiling((n - 1) / 2)
    ) %>%
      dplyr::mutate(
        lower = pmax(
          min_lower,
          fitted - pmax(min_radius, detection_multiplier * roll_iqr)
        ),
        upper = fitted + pmax(min_radius, detection_multiplier * roll_iqr),
        replacement = dplyr::case_when(
          (y < lower) ~ as_type(fitted - replacement_multiplier * roll_iqr),
          (y > upper) ~ as_type(fitted + replacement_multiplier * roll_iqr),
          TRUE ~ y
        )
      ) %>%
      dplyr::select(lower, upper, replacement) %>%
      tibble::as_tibble()
  }
  z = as_epi_df(tibble::tibble(
    geo_value = 0,
    time_value = seq_along(deaths),
    y = deaths
  ))
  epidf <- epi_slide(
    z,
    fitted = median(y),
    before = floor((n - 1) / 2),
    after = ceiling((n - 1) / 2)
  )
  s <- dplyr::mutate(epidf, resid = y - fitted) # make new column, resid
  z2 <- roll_iqr(s, detection_multiplier = 5)$replacement
  return(z2)
}

####################

# make_delay_distr <- function(Mean, Sd, d) {
#   params <- get_gamma_params(Mean=Mean, Sd=Sd)
#   shape <- params[1]; rate <- params[2]
#   delay_distr <- ddgamma(0:d, shape, rate)
#   return(delay_distr)
# }

smooth_counts <- function(x, w, centered = TRUE) {
  nm <- names(x)
  x <- unname(x)
  if (centered) {
    l <- floor((w - 1) / 2)
    r <- floor(w / 2)
    z <- c(
      tail(drop(lower.tri(matrix(NA, w, w)) %*% head(x, w)) / 0:(w - 1), l),
      roll_mean(x, n = w),
      head(drop(upper.tri(matrix(NA, w, w)) %*% tail(x, w)) / (w - 1):0, r)
    )
  } else {
    z <- roll_meanr(x, n = w, fill = mean(head(x, w)))
  }
  if (!is_null(nm)) {
    names(z) <- nm
  }
  z
}

# smooth_counts_old <- function(data, w, centered = TRUE) {
#   n <- length(data)
#   if (w <= 1) {
#     return(data)
#   }
#   if (centered == FALSE) {
#     smoothed <- sapply(1:n, function(i) {
#       # Will be flat for the first 7 days
#       first <- ifelse(i < w, 1, i - w + 1)
#       last <- ifelse(i < w, w, i)
#       mean(data[first:last])
#     })
#   } else {
#     smoothed <- sapply(1:n, function(i) {
#       # Later, should account for boundary effects
#       n_into_past <- floor((w - 1) / 2)
#       n_into_future <- floor(w / 2)
#       first <- ifelse(i - n_into_past < 1, 1, i - n_into_past)
#       last <- ifelse(i + n_into_future > n, n, i + n_into_future)
#       mean(data[first:last])
#     })
#   }
#   if (!is.null(names(data))) {
#     names(smoothed) <- names(data)
#   }
#   return(smoothed)
# }

clamp <- function(x, ll = 0, uu = 1) {
  pmax(pmin(x, uu), ll)
}
