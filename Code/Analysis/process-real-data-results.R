library(tidyverse)
source(here::here("Code", "funs.R"))
source(here::here("Code", "estimate.R"))
source(here::here("Code", "helper_functions.R"))
source(here::here("Code", "Analysis", "collection-funs.R"))

res <- read_rds(here::here("Data", "batch-results", "real-data-res.rds"))
meta <- res$meta
region <- meta$region[1]

# Duplicate some processing to get the correct retro time_value
data_finalized <- read_rds(here::here("Data", "raw-data", "finalized.rds")) |>
  filter(geo_value == region)
nchs_daily <- data_finalized |>
  select(-hosps) |>
  filter(wday(time_value) == 7) |>
  rowwise() |>
  mutate(tvmod = list(6:0)) |>
  ungroup() |>
  unnest_longer(tvmod) |>
  mutate(time_value = ymd(time_value - tvmod), tvmod = NULL)
time_value <- inner_join(
  data_finalized,
  nchs_daily,
  by = c("geo_value", "time_value")
) |>
  arrange(time_value) |>
  pull(time_value)


# retro -----------------------------------------------------------------------
filter_intersect_tv <- function(tib) {
  n_levels <- n_distinct(tib$estimator)
  tib |>
    group_by(time_value) |>
    filter(n_distinct(estimator) == n_levels) |>
    ungroup()
}

tuned_retro <- function(res, which_tune = c("min", "1se"), benchmark = TRUE) {
  which_tune <- match.arg(which_tune)
  idx <- paste("i", which_tune, sep = ".")
  idx <- res$cvstats[[idx]]
  if (benchmark) {
    fits <- res$full_fit[, idx] |>
      enframe("time_value", "est") |>
      mutate(time_value = as.Date(time_value))
  } else {
    fits <- tibble(time_value = time_value, est = res$full_fit[, idx])
  }
  fits
}
retro <- c(res[[1]], res[2:4] |> set_names(paste0("tf", 0:2)))
retro_min <- map2(retro, c(TRUE, TRUE, rep(FALSE, 3)), \(x, y) {
  tuned_retro(x, "min", y)
}) |>
  list_rbind(names_to = "estimator") |>
  filter_intersect_tv()

retro_1se <- map2(retro, c(TRUE, TRUE, rep(FALSE, 3)), \(x, y) {
  tuned_retro(x, "1se", y)
}) |>
  list_rbind(names_to = "estimator") |>
  filter_intersect_tv()

# ggplot(retro_min, aes(time_value, est, color = estimator)) +
#   geom_line() +
#   theme_bw()

# real time -------------------------------------------------------------------
rt <- c(res[5], res[6:8] |> set_names(paste0("tf", 0:2)))
rt_curves <- function(rt_list, cv = c("min", "1se"), fv = c("min", "1se")) {
  lagged <- benchcurve_rt(rt_list[[1]], 1, fv, "lagged")
  conv <- benchcurve_rt(rt_list[[1]], 1, fv, "conv")
  tf <- map(rt[2:4], ~ tune_tf_order(.x, cv, fv) |> select(-cverr))
  c(list(lagged = lagged, conv = conv), tf) |>
    list_rbind(names_to = "estimator")
}
rt_min_min <- rt_curves(rt, "min", "min")
rt_min_1se <- rt_curves(rt, "min", "1se")
rt_1se_min <- rt_curves(rt, "1se", "min")
rt_1se_1se <- rt_curves(rt, "1se", "1se")

real_data_pa_curves <- enlist(
  retro_min,
  retro_1se,
  rt_min_min,
  rt_min_1se,
  rt_1se_min,
  rt_1se_1se
)

write_rds(
  real_data_pa_curves,
  here::here("Data", "batch-results", "real-data-curves.rds")
)
