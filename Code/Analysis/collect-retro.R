library(batchtools)
library(tidyverse)
loadRegistry("severity-estimation-retro-v2")
burn_len <- 60L
hfrs <- read_rds(here::here("Data", "raw-data", "daily-hfrs.rds"))

source(here::here("Code", "Analysis", "collection-funs.R"))

allids <- getJobTable() |>
  select(job.id, repl, problem, prob.pars, algorithm, algo.pars) |>
  unwrap()
ids_tf <- allids |> filter(algorithm == "trendfilter")
ids_b <- allids |> filter(algorithm == "benchmarks")

# all tuned curves ------------------------------------------------------------
tfcurves <- reduceResultsDataTable(
  ids = ids_tf,
  fun = tfcurve_retro,
  which_lambda = "min"
) |>
  left_join(ids_tf, by = "job.id")
convcurves <- reduceResultsDataTable(
  ids_b,
  fun = benchcurve_retro,
  which_w = "1se",
  estimator = "conv"
) |>
  left_join(ids_b, by = "job.id")
lagcurves <- reduceResultsDataTable(
  ids_b,
  fun = benchcurve_retro,
  which_w = "1se",
  estimator = "lagged"
) |>
  left_join(ids_b, by = "job.id")

allres <- bind_rows(
  trendfilter = tfcurves,
  conv = convcurves,
  lagged = lagcurves,
  .id = "estimator"
)

cverrs_b <- reduceResultsDataTable(ids_b, fun = function(job, res) {
  tibble(
    conv = res$conv$cvstats$cvm[res$conv$cvstats$i.min],
    lagged = res$lagged$cvstats$cvm[res$lagged$cvstats$i.min]
  )
}) |>
  unnest(result) |>
  pivot_longer(conv:lagged, names_to = "estimator", values_to = "cverr")
cverrs_tf <- reduceResultsDataTable(ids_tf, fun = function(job, res) {
  tibble(cverr = res$cvstats$cvm[res$cvstats$i.min])
}) |>
  unnest(result) |>
  mutate(estimator = "trendfilter")
cverrs <- bind_rows(cverrs_b, cverrs_tf) |>
  left_join(allids, by = join_by(job.id)) |>
  select(estimator, cverr, repl, geo_value = region, mod_idx, order)

tuned_maes <- allres |>
  unnest(result) |>
  select(
    job.id,
    time_value,
    est,
    repl,
    geo_value = region,
    mod_idx,
    order,
    estimator
  ) |>
  arrange(job.id, estimator, order, time_value) |>
  slice(
    -(1:burn_len),
    -((n() - burn_len):n()),
    .by = c(job.id, estimator, order)
  ) |> # burn-in/-out
  left_join(hfrs, by = join_by(time_value, geo_value)) |>
  group_by(repl, geo_value, mod_idx) |>
  group_modify(
    ~ {
      common_time_values <- Reduce(
        intersect,
        split(.x$time_value, paste0(.x$estimator, .x$order))
      )
      .x |>
        filter(time_value %in% common_time_values) |>
        arrange(estimator, order, time_value)
    }
  ) |>
  ungroup() |>
  group_by(repl, geo_value, mod_idx, order, estimator) |>
  summarise(
    mae = mean(abs(replace_inf(est) - hfr), na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(cverrs)


saveRDS(
  tuned_maes |> filter(mod_idx == 1),
  here::here("Data", "batch-results", "retro-maes-ws.rds")
)
saveRDS(
  tuned_maes |> filter(repl == 1),
  here::here("Data", "batch-results", "retro-maes-mis.rds")
)
saveRDS(
  allres |> filter(mod_idx == 1),
  here::here("Data", "batch-results", "retro-example-curves.rds")
)

# oracle curves ---------------------------------------------------------------
tfcurves <- reduceResultsDataTable(
  ids = ids_tf,
  fun = tfcurve_retro,
  which_lambda = "oracle"
) |>
  left_join(ids_tf, by = "job.id")
convcurves <- reduceResultsDataTable(
  ids_b,
  fun = benchcurve_retro,
  which_w = "oracle",
  estimator = "conv"
) |>
  left_join(ids_b, by = "job.id")
lagcurves <- reduceResultsDataTable(
  ids_b,
  fun = benchcurve_retro,
  which_w = "oracle",
  estimator = "lagged"
) |>
  left_join(ids_b, by = "job.id")

oracle_curves <- bind_rows(
  trendfilter = tfcurves,
  conv = convcurves,
  lagged = lagcurves,
  .id = "estimator"
) |>
  filter(mod_idx == 1)

oracle_maes <- oracle_curves |>
  unnest(result) |>
  select(
    job.id,
    time_value,
    est,
    repl,
    geo_value = region,
    mod_idx,
    order,
    estimator
  ) |>
  arrange(job.id, estimator, order, time_value) |>
  slice(
    -(1:burn_len),
    -((n() - burn_len):n()),
    .by = c(job.id, estimator, order)
  ) |> # burn-in/-out
  left_join(hfrs, by = join_by(time_value, geo_value)) |>
  group_by(repl, geo_value, mod_idx) |>
  group_modify(
    ~ {
      common_time_values <- Reduce(
        intersect,
        split(.x$time_value, paste0(.x$estimator, .x$order))
      )
      .x |>
        filter(time_value %in% common_time_values) |>
        arrange(estimator, order, time_value)
    }
  ) |>
  ungroup() |>
  group_by(repl, geo_value, mod_idx, order, estimator) |>
  summarise(
    mae = mean(abs(replace_inf(est) - hfr), na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(cverrs)

saveRDS(
  oracle_maes,
  here::here("Data", "batch-results", "retro-maes-ws-oracle.rds")
)
