library(batchtools)
source(here::here("Code", "Analysis", "collection-funs.R"))

loadRegistry("severity-estimation-rt-v2")

gt_hfrs <- read_rds(here::here("Data", "raw-data", "daily-hfrs.rds")) |>
  select(geo_value, time_value, hfr)


allids <- getJobTable() |>
  select(job.id, repl, problem, prob.pars, algorithm, algo.pars) |>
  unwrap()
wsids <- allids |> filter(mod_idx == 1)
msids <- allids |> filter(repl == 1)
ids_tf <- wsids |> filter(algorithm == "trendfilter")
ids_b <- wsids |> filter(algorithm == "benchmarks")

tferrs <- reduceResultsDataTable(
  ids_tf,
  grab_rt_maes,
  cvopt = "min",
  fvopt = "1se"
) |>
  unwrap() |>
  mutate(estimator = "trendfilter", fv = "1se", cv = "min")

converrs <- reduceResultsDataTable(
  ids_b,
  fun = grab_rt_maes,
  fvopt = "1se",
  estimator = "conv"
) |>
  unwrap() |>
  mutate(estimator = "conv")
lagerrs <- reduceResultsDataTable(
  ids_b,
  fun = grab_rt_maes,
  fvopt = "1se",
  estimator = "lagged"
) |>
  unwrap() |>
  mutate(estimator = "lagged")

tf_curves_to_tune_order <- reduceResultsList(
  ids_tf,
  tune_tf_order,
  cvopt = "min",
  fvopt = "1se"
)
tf_to_tune_dt <- getJobTable(ids_tf) |>
  select(repl, algo.pars, prob.pars) |>
  unwrap()
tf_to_tune_dt$fits <- tf_curves_to_tune_order

tf_tuned_maes <- tf_to_tune_dt |>
  group_by(region, repl) |>
  group_modify(
    ~ {
      local_hfr <- gt_hfrs |> filter(geo_value == .y$region[1])
      z <- .x |> unnest(fits)
      tune_idx <- summarise(
        z,
        best_order = min_order(order, est, cverr),
        .by = c(time_value)
      )
      z |>
        select(order, time_value, est) |>
        mutate(order = factor(order, levels = 0:2)) |>
        pivot_wider(
          names_from = order,
          values_from = est,
          names_expand = TRUE
        ) |>
        left_join(tune_idx, by = "time_value") |>
        transmute(
          time_value = time_value,
          est = case_when(
            best_order == 0 ~ `0`,
            best_order == 1 ~ `1`,
            best_order == 2 ~ `2`,
            TRUE ~ `0`
          )
        ) |>
        left_join(local_hfr, by = "time_value") |>
        arrange(time_value) |>
        fill(hfr, .direction = "downup") |>
        slice(-(1:3)) |>
        summarise(mae = mean(abs(est - hfr), na.rm = TRUE))
    }
  ) |>
  ungroup() |>
  rename(geo_value = region) |>
  mutate(estimator = "tf_tuned_order")

tf_tuned_maes <- tf_tuned_maes |> mutate(fv = "1se", cv = "min")
allres <- bind_rows(tferrs, converrs, lagerrs, tf_tuned_maes)
saveRDS(allres, here::here("Data", "batch-results", "rt-maes-ws.rds"))


# Curves ------------------------------------------------------------------

tfcurves <- tf_to_tune_dt |>
  group_by(region, order) |>
  # slice_head(n = 1) |>
  ungroup() |>
  unnest(fits) |>
  select(order, region, repl, time_value, est) |>
  mutate(estimator = "trendfilter")
tf_tuned_curves <- tf_to_tune_dt |>
  # group_by(region, order) |>
  # slice_head(n = 1) |>
  # ungroup() |>
  group_by(region, repl) |>
  group_modify(
    ~ {
      z <- .x |> unnest(fits)
      tune_idx <- summarise(
        z,
        best_order = min_order(order, est, cverr),
        .by = c(time_value)
      )
      z |>
        select(order, time_value, est) |>
        mutate(order = factor(order, levels = 0:2)) |>
        pivot_wider(
          names_from = order,
          values_from = est,
          names_expand = TRUE
        ) |>
        left_join(tune_idx, by = "time_value") |>
        transmute(
          time_value = time_value,
          est = case_when(
            best_order == 0 ~ `0`,
            best_order == 1 ~ `1`,
            best_order == 2 ~ `2`,
            TRUE ~ `0`
          )
        )
    }
  ) |>
  ungroup() |>
  mutate(estimator = "tf_tuned_order")


convcurves <- ids_b |> # slice_head(n = 1, by = region) |>
  reduceResultsDataTable(fun = benchcurve_rt) |>
  left_join(ids_b |> select(job.id, region, repl), by = "job.id") |>
  mutate(estimator = "conv") |>
  select(-job.id) |>
  unnest(result)
laggedcurves <- ids_b |> # slice_head(n = 1, by = region) |>
  reduceResultsDataTable(fun = benchcurve_rt, estimator = "lagged") |>
  left_join(ids_b |> select(job.id, region, repl), by = "job.id") |>
  mutate(estimator = "lagged") |>
  select(-job.id) |>
  unnest(result)
allcurves <- bind_rows(tfcurves, tf_tuned_curves, convcurves, laggedcurves)

saveRDS(allcurves, here::here("Data", "batch-results", "rt-example-curves.rds"))


# Misspecified ------------------------------------------------------------

# One VT did not complete
bad_ids <- findErrors()$job.id
ids_tf <- msids |> filter(algorithm == "trendfilter", job.id %nin% bad_ids)
ids_b <- msids |> filter(algorithm == "benchmarks")

tferrs <- reduceResultsDataTable(
  ids_tf,
  grab_rt_maes,
  cvopt = "min",
  fvopt = "1se"
) |>
  unwrap() |>
  mutate(estimator = "trendfilter", fv = "1se", cv = "min") |>
  left_join(ids_tf |> select(job.id, mod_idx))
converrs <- reduceResultsDataTable(
  ids_b,
  fun = grab_rt_maes,
  fvopt = "1se",
  estimator = "conv"
) |>
  unwrap() |>
  mutate(estimator = "conv") |>
  left_join(ids_b |> select(job.id, mod_idx))
lagerrs <- reduceResultsDataTable(
  ids_b,
  fun = grab_rt_maes,
  fvopt = "1se",
  estimator = "lagged"
) |>
  unwrap() |>
  mutate(estimator = "lagged") |>
  left_join(ids_b |> select(job.id, mod_idx))

tf_curves_to_tune_order <- reduceResultsList(
  ids_tf,
  tune_tf_order,
  cvopt = "min",
  fvopt = "1se"
)
tf_to_tune_dt <- getJobTable(ids_tf) |>
  select(repl, algo.pars, prob.pars) |>
  unwrap()
tf_to_tune_dt$fits <- tf_curves_to_tune_order

tf_tuned_maes <- tf_to_tune_dt |>
  group_by(region, mod_idx) |>
  group_modify(
    ~ {
      local_hfr <- gt_hfrs |> filter(geo_value == .y$region[1])
      z <- .x |> unnest(fits)
      tune_idx <- summarise(
        z,
        best_order = min_order(order, est, cverr),
        .by = c(time_value)
      )
      z |>
        select(order, time_value, est) |>
        mutate(order = factor(order, levels = 0:2)) |>
        pivot_wider(
          names_from = order,
          values_from = est,
          names_expand = TRUE
        ) |>
        left_join(tune_idx, by = "time_value") |>
        transmute(
          time_value = time_value,
          est = case_when(
            best_order == 0 ~ `0`,
            best_order == 1 ~ `1`,
            best_order == 2 ~ `2`,
            TRUE ~ `0`
          )
        ) |>
        left_join(local_hfr, by = "time_value") |>
        arrange(time_value) |>
        fill(hfr, .direction = "downup") |>
        slice(-(1:3)) |>
        summarise(mae = mean(abs(est - hfr), na.rm = TRUE))
    }
  ) |>
  ungroup() |>
  rename(geo_value = region) |>
  mutate(estimator = "tf_tuned_order")

tf_tuned_maes <- tf_tuned_maes |> mutate(fv = "1se", cv = "min")


allres <- bind_rows(tferrs, lagerrs, converrs)
allres <- bind_rows(allres, tf_tuned_maes)
saveRDS(allres, here::here("Data", "batch-results", "rt-maes-mis.rds"))
