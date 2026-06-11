library(batchtools)
library(tidyverse)
library(data.table)
source(here::here("Code", "helper_functions.R"))
source(here::here("Code", "estimate.R"))
source(here::here("Code", "funs.R"))
source(here::here("Code", "batch-funs.R"))

raw_data <- read_rds(here::here("Data", "raw-data", "finalized.rds")) |>
  select(-deaths)
gt_hfrs <- read_rds(here::here("Data", "raw-data", "daily-hfrs.rds"))
meta_data <- read_rds(here::here("Data", "raw-data", "all-meta-data.rds"))

regions <- sort(meta_data$geo_value) # Includes US

## Global parameters
d <- 60L
ws <- c(1, 7, 15, 21, 29)
ord <- 0:2

## Delay mean modifier for misspecified retrospective case.
delay_mean_modifier <- c(0, 0, 0, -3, -2, -1, 1, 2, 3)
delay_sds <- c(.9, 1, .8, rep(.9, 6))

data_list <- enlist(
  raw_data,
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
  benchmarks = data.frame(mod_idx = 2:9)
)
algo_design_ws <- list(
  trendfilter = data.frame(mod_idx = 1, order = c(0, 1, 2)),
  benchmarks = data.frame(mod_idx = 1)
)


# Set up the experimental design
makeExperimentRegistry(
  "severity-estimation-retro-v2",
  source = c(
    here::here("Code", "helper_functions.R"),
    here::here("Code", "estimate.R"),
    here::here("Code", "funs.R"),
    here::here("Code", "batch-funs.R")
  ),
  packages = c("tidyverse")
)
addProblem(
  "regional_hfr",
  fun = create_regional_problem_retro,
  data = data_list,
  seed = 123
)
addAlgorithm("trendfilter", fun = trendfilter_retro_algo)
addAlgorithm("benchmarks", fun = benchmark_retro_algos)
addExperiments(prob_design, algo_design_misp)
addExperiments(prob_design, algo_design_ws, repls = 10L)

ids <- getJobPars() |>
  unwrap() |>
  mutate(chunk = as.numeric(as.factor(region))) |>
  select(job.id, chunk)
# submitJobs(ids, resources = list(walltime = "12:0:0", memory = "4gb", ncpus = 1))
