# Real time functions ---------------------------------------------------------
create_regional_problem_rt <- function(data, region, time_var = FALSE, ...) {
  local_meta <- data$meta_data |> filter(geo_value == region)
  local_data <- data$data_archive |>
    filter(geo_value == region)
  tvar_delay <- NULL
  if (is_null(time_var)) {
    time_var <- FALSE
  }
  time_var <- time_var & !is_null(data$variant_props)
  delay_mean <- local_meta |> pull(rt_lag_1se)
  delay_distn <- make_delay_distr(delay_mean, 0.9 * delay_mean, data$d)
  dispersion <- max(local_meta |> pull(dispersion), 1)
  sim_partial <- function(x, ndates = 7) {
    ll <- process_tibble_to_named_vctrs(x)
    if (time_var) {
      tvar_delay <- x$rt_delay
      tvar_delay <- do.call(rbind, tvar_delay)
    }
    if (dispersion > 1) {
      res <- sim_bb_deaths(
        ll$death_dates,
        ll$hh,
        ll$hfr,
        delay_distn,
        dispersion,
        tvar_delay = tvar_delay
      )
    } else {
      res <- sim_pb_deaths(
        ll$death_dates,
        ll$hh,
        ll$hfr,
        delay_distn,
        tvar_delay = tvar_delay
      )
    }
    res <- enframe(res, "time_value", "sim_deaths") |>
      mutate(time_value = ymd(time_value), geo_value = region) |>
      arrange(time_value) |>
      slice_tail(n = ndates)
  }
  first_version <- min(local_data$DT$version)
  res <- local_data |>
    epix_slide(.f = function(x, g, v) {
      if (time_var) {
        x <- x |>
          left_join(data$variant_props, by = join_by(geo_value, time_value))
      }
      x |>
        inner_join(
          data$gt_hfrs |> filter(geo_value == region),
          by = join_by(geo_value, time_value)
        ) |>
        arrange(time_value) |>
        sim_partial(ndates = ifelse(v > first_version, 7, Inf))
    }) |>
    as.data.frame() |>
    arrange(version, time_value) |>
    as_epi_archive()
  res$geo_type <- "nation"
  res
}

benchmark_rt_algos <- function(
  job,
  data,
  instance,
  mod_idx = NULL,
  spec = NULL,
  ...
) {
  region <- head(instance$DT$geo_value, 1)
  local_meta <- data$meta_data |> filter(geo_value == region)
  local_data <- data$data_archive |>
    filter(geo_value == region) |>
    epix_merge(instance, sync = "truncate")
  local_data$DT <- local_data$DT |>
    select(-deaths) |>
    rename(deaths = sim_deaths)
  oracle_delay_mean <- local_meta |> pull(rt_lag_1se)
  if (is_null(mod_idx)) {
    delay_mean <- oracle_delay_mean
    delay_sd <- delay_mean * 0.9
  } else {
    delay_mean <- max(oracle_delay_mean + data$delay_mean_modifier[mod_idx], 1)
    delay_sd <- delay_mean * data$delay_sds[mod_idx]
  }
  delay_distn <- make_delay_distr(delay_mean, delay_sd, data$d)
  est_dates <- sort(unique(local_data$DT$version))
  lagged <- epix_slide(
    local_data,
    .versions = est_dates,
    .f = function(x, g, v) {
      tvar_delay <- NULL
      tvar_lag <- NULL
      if (!is_null(spec) && !is_null(data$variant_props) && spec == "ws") {
        tvar_delay <- x |>
          left_join(data$variant_props, by = join_by(geo_value, time_value)) |>
          pull(rt_delay)
        tvar_delay <- do.call(rbind, tvar_delay)
        tvar_lag <- pmax(pmin(round(drop(tvar_delay %*% 0:data$d)), 35), 5)
      }
      ll <- process_tibble_to_named_vctrs(x)
      list(tune_lagged_hfrs_rt(
        hosps = ll$hh,
        deaths = ll$dd,
        t = v,
        delay_distr = delay_distn,
        ws = data$ws,
        lag = delay_mean,
        tvar_lag = tvar_lag,
        tvar_delay = tvar_delay
      ))
    }
  )
  conv <- epix_slide(
    local_data,
    .versions = est_dates,
    .f = function(x, g, v) {
      tvar_delay <- NULL
      if (!is_null(spec) && !is_null(data$variant_props) && spec == "ws") {
        tvar_delay <- x |>
          left_join(data$variant_props, by = join_by(geo_value, time_value)) |>
          pull(rt_delay)
        tvar_delay <- do.call(rbind, tvar_delay)
      }
      ll <- process_tibble_to_named_vctrs(x)
      list(tune_conv_hfrs_rt(
        hosps = ll$hh,
        deaths = ll$dd,
        t = v,
        delay_distr = delay_distn,
        ws = data$ws,
        tvar_delay = tvar_delay
      ))
    }
  )
  enlist(lagged, conv)
}

trendfilter_rt_algo <- function(
  job,
  data,
  instance,
  mod_idx = NULL,
  order = 0,
  spec = NULL,
  ...
) {
  gammas <- geom_spacing(1, 1e5, 20L)
  region <- head(instance$DT$geo_value, 1)
  local_meta <- data$meta_data |> filter(geo_value == region)
  local_data <- data$data_archive |>
    filter(geo_value == region) |>
    epix_merge(instance, sync = "truncate")
  local_data$DT <- local_data$DT |>
    select(-deaths) |>
    rename(deaths = sim_deaths)
  oracle_delay_mean <- local_meta |> pull(rt_lag_1se)
  if (is_null(mod_idx)) {
    delay_mean <- oracle_delay_mean
    delay_sd <- 0.9 * oracle_delay_mean
  } else {
    delay_mean <- max(oracle_delay_mean + data$delay_mean_modifier[mod_idx], 1)
    delay_sd <- delay_mean * data$delay_sds[mod_idx]
  }
  delay_distn <- make_delay_distr(delay_mean, delay_sd, data$d)
  est_dates <- sort(unique(local_data$DT$version))
  res <- epix_slide(
    local_data,
    .versions = est_dates,
    .f = function(x, g, v) {
      tvar_delay <- NULL
      if (!is_null(spec) && !is_null(data$variant_props) && spec == "ws") {
        tvar_delay <- x |>
          left_join(data$variant_props, by = join_by(geo_value, time_value)) |>
          pull(rt_delay)
        tvar_delay <- do.call(rbind, tvar_delay)
      }
      ll <- process_tibble_to_named_vctrs(x)
      list(tune_trendfilter_rt(
        ll$hh,
        ll$dd,
        gammas,
        delay_distn,
        order = order,
        tvar_delay = tvar_delay
      ))
    }
  )
  res
}

# Retrospective functions -----------------------------------------------------
create_regional_problem_retro <- function(data, region, time_var = FALSE, ...) {
  local_meta <- data$meta_data |> filter(geo_value == region)
  local_data <- data$raw_data |>
    filter(geo_value == region) |>
    left_join(data$gt_hfrs, by = join_by(geo_value, time_value))
  tvar_delay <- NULL
  if (is_null(time_var)) {
    time_var <- FALSE
  }
  time_var <- time_var & !is_null(data$variant_props)
  if (time_var) {
    tvar_delay <- local_data |>
      left_join(data$variant_props, by = join_by(time_value, geo_value)) |>
      pull(retro_delay)
    tvar_delay <- do.call(rbind, tvar_delay)
  }
  delay_mean <- local_meta |> pull(retro_lag_argmax)
  delay_distn <- make_delay_distr(delay_mean, 0.9 * delay_mean, data$d)
  hh <- with(local_data, set_names(hosps, time_value))
  hfr <- with(local_data, set_names(hfr, time_value))
  death_dates <- local_data$time_value[-c(1:data$d)]
  sim_pb_deaths(death_dates, hh, hfr, delay_distn, tvar_delay = tvar_delay) |>
    enframe("time_value", "sim_deaths") |>
    mutate(time_value = ymd(time_value), geo_value = region) |>
    arrange(time_value)
}

benchmark_retro_algos <- function(
  job,
  data,
  instance,
  mod_idx = NULL,
  spec = NULL,
  ...
) {
  region <- head(instance$geo_value, 1)
  local_meta <- data$meta_data |> filter(geo_value == region)
  local_data <- data$raw_data |>
    filter(geo_value == region) |>
    left_join(instance, by = join_by(geo_value, time_value)) |>
    rename(deaths = sim_deaths)
  oracle_delay_mean <- local_meta |> pull(retro_lag_argmax)
  if (is_null(mod_idx)) {
    delay_mean <- oracle_delay_mean
    delay_sd <- delay_mean * 0.9
  } else {
    delay_mean <- max(oracle_delay_mean + data$delay_mean_modifier[mod_idx], 1)
    delay_sd <- delay_mean * data$delay_sds[mod_idx]
  }
  delay_distn <- make_delay_distr(delay_mean, delay_sd, data$d)
  ll <- process_tibble_to_named_vctrs(local_data)
  tvar_delay <- NULL
  tvar_lag <- NULL
  if (!is_null(spec) && !is_null(data$variant_props) && spec == "ws") {
    tvar_delay <- local_data |>
      left_join(data$variant_props, by = join_by(geo_value, time_value)) |>
      pull(retro_delay)
    tvar_delay <- do.call(rbind, tvar_delay)
    tvar_lag <- pmax(pmin(round(drop(tvar_delay %*% 0:data$d)), 35), 5)
  }
  lagged <- tune_lagged_hfrs_retro(
    hosps = ll$hh,
    deaths = ll$dd,
    delay_distr = delay_distn,
    ws = data$ws,
    lag = delay_mean,
    tvar_lag = tvar_lag,
    tvar_delay = tvar_delay
  )
  conv <- tune_conv_hfrs_retro(
    ll$hh,
    ll$dd,
    delay_distn,
    data$ws,
    tvar_delay = tvar_delay
  )
  enlist(lagged, conv)
}

trendfilter_retro_algo <- function(
  job,
  data,
  instance,
  mod_idx = NULL,
  order = 0,
  spec = NULL,
  ...
) {
  region <- head(instance$geo_value, 1)
  local_meta <- data$meta_data |> filter(geo_value == region)
  local_data <- data$raw_data |>
    filter(geo_value == region) |>
    left_join(instance, by = join_by(geo_value, time_value)) |>
    rename(deaths = sim_deaths)
  oracle_delay_mean <- local_meta |> pull(retro_lag_argmax)
  if (is_null(mod_idx)) {
    delay_mean <- oracle_delay_mean
    delay_sd <- 0.9 * oracle_delay_mean
  } else {
    delay_mean <- max(oracle_delay_mean + data$delay_mean_modifier[mod_idx], 1)
    delay_sd <- delay_mean * data$delay_sds[mod_idx]
  }
  tvar_delay <- NULL
  delay_distn <- make_delay_distr(delay_mean, delay_sd, data$d)
  if (!is_null(spec) && !is_null(data$variant_props) && spec == "ws") {
    tvar_delay <- local_data |>
      left_join(data$variant_props, by = join_by(geo_value, time_value)) |>
      pull(retro_delay)
    tvar_delay <- do.call(rbind, tvar_delay)
  }
  ll <- process_tibble_to_named_vctrs(local_data)
  cv_lambda_tf(
    ll$hh,
    ll$dd,
    delay_distn,
    order = order,
    tvar_delay = tvar_delay
  )
}
