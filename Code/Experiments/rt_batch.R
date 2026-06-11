library(batchtools)
library(tidyverse)
library(epiprocess)
library(data.table)
source(here::here("Code", "helper_functions.R"))
source(here::here("Code", "estimate.R"))
source(here::here("Code", "funs.R"))
source(here::here("Code", "batch-funs.R"))

data_archive <- read_rds(here::here("Data", "raw-data", "full-archive.rds")) |>
  as_epi_archive()
gt_hfrs <- read_rds(here::here("Data", "raw-data", "daily-hfrs.rds"))
meta_data <- read_rds(here::here("Data", "raw-data", "all-meta-data.rds"))

regions <- sort(meta_data$geo_value)

## Global parameters
d <- 60L
ws <- c(1, 7, 15, 21, 29)
ord <- 0:2

## First condition is correctly specified
delay_mean_modifier <- c(0, 0, 0, -5, -3, -1, 1, 3, 5)
delay_sds <- c(.9, 1, .8, rep(.9, 6))

data_list <- enlist(
  data_archive,
  gt_hfrs,
  meta_data,
  delay_mean_modifier,
  delay_sds,
  ws,
  d
)

prob_design <- list(regional_hfr = data.frame(region = regions))
algo_design_misp <- list(
  trendfilter = CJ(mod_idx = 2:9, order = c(0, 1, 2)),
  benchmarks = data.frame(mod_idx = seq_along(delay_mean_modifier)[-1])
)
algo_design_ws <- list(
  trendfilter = data.frame(mod_idx = 1, order = c(0, 1, 2)),
  benchmarks = data.frame(mod_idx = 1)
)

# Set up the experimental design
makeExperimentRegistry(
  "severity-estimation-rt-v2",
  source = c(
    here::here("Code", "helper_functions.R"),
    here::here("Code", "estimate.R"),
    here::here("Code", "funs.R"),
    here::here("Code", "batch-funs.R")
  ),
  packages = c("tidyverse", "epiprocess")
)
addProblem(
  "regional_hfr",
  fun = create_regional_problem_rt,
  data = data_list,
  seed = 123
)
addAlgorithm("trendfilter", fun = trendfilter_rt_algo)
addAlgorithm("benchmarks", fun = benchmark_rt_algos)
addExperiments(prob_design, algo_design_misp)
addExperiments(prob_design, algo_design_ws, repls = 10L)

ids <- getJobTable() |>
  unwrap() |>
  select(job.id, repl, problem, algorithm, region, mod_idx, order) |>
  arrange(region, mod_idx, repl) |>
  mutate(chunk = rep(seq_len(n() / 8), each = 8)) |>
  arrange(job.id)
# submitJobs(ids, resources = list(walltime = "48:0:0", memory = "4gb", ncpus = 1))
