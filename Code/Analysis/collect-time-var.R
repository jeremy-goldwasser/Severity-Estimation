library(batchtools)
source(here::here("Code", "Analysis", "collection-funs.R"))

loadRegistry("severity-estimation-tvar-delay")

allids <- getJobTable() |>
  select(job.id, repl, problem, prob.pars, algorithm, algo.pars) |>
  unwrap()
ids_tf_rt <- allids |> filter(algorithm == "tf-rt")
ids_b_rt <- allids |> filter(algorithm == "bench-rt")
ids_tf_retro <- allids |> filter(algorithm == "tf-retro")
ids_b_retro <- allids |> filter(algorithm == "bench-retro")

ids_tf_rt <- ids_tf_rt |>
  mutate(
    res = reduceResultsList(
      ids = job.id,
      fun = tune_tf_order,
      cvopt = "min",
      fvopt = "1se"
    ),
    estimator = "trendfilter"
  )
ids_tf_retro <- ids_tf_retro |>
  mutate(
    res = reduceResultsList(
      ids = job.id,
      fun = tfcurve_retro,
      which_lambda = "min"
    ),
    estimator = "trendfilter"
  )
ids_b_rt <- ids_b_rt |>
  mutate(
    conv = reduceResultsList(ids = job.id, fun = benchcurve_rt),
    lagged = reduceResultsList(
      ids = job.id,
      fun = benchcurve_rt,
      estimator = "lagged"
    )
  ) |>
  pivot_longer(c(conv, lagged), names_to = "estimator", values_to = "res")
ids_b_retro <- ids_b_retro |>
  mutate(
    conv = reduceResultsList(
      ids = job.id,
      fun = benchcurve_retro,
      which_w = "1se",
      estimator = "conv"
    ),
    lagged = reduceResultsList(
      ids = job.id,
      fun = benchcurve_retro,
      which_w = "1se",
      estimator = "lagged"
    )
  ) |>
  pivot_longer(c(conv, lagged), names_to = "estimator", values_to = "res")


res <- bind_rows(ids_tf_rt, ids_b_rt, ids_tf_retro, ids_b_retro)

write_rds(res, here::here("Data", "batch-results", "time-var-delay.rds"))
