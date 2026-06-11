library(batchtools)
library(tidyverse)
library(epiprocess)
library(data.table)

source(here::here("Code", "funs.R"))
source(here::here("Code", "helper_functions.R"))
source(here::here("Code", "estimate.R"))
source(here::here("Code", "batch-funs.R"))


# Global parameters -------------------------------------------------------

d <- 60L
ws <- c(1, 7, 15, 21, 29)
ord <- 0:2
region <- "pa"

# Grab data ---------------------------------------------------------------
data_finalized <- read_rds(here::here("Data", "raw-data", "finalized.rds")) |>
  filter(geo_value == region)
data_archive <- read_rds(here::here("Data", "raw-data", "full-archive.rds")) |>
  as.data.frame() |>
  filter(geo_value == region) |>
  as_epi_archive()
meta_data <- read_rds(here::here("Data", "raw-data", "all-meta-data.rds"))

# Fix retro weekly death reporting
set.seed(12345)
nchs_daily <- data_finalized |>
  select(-hosps) |>
  filter(wday(time_value) == 7) |>
  mutate(deaths = deaths * 7) |>
  rowwise() |>
  mutate(tvmod = list(6:0), deaths = list(rmultinom1(deaths, rep(1, 7)))) |>
  ungroup() |>
  unnest_longer(c(tvmod, deaths)) |>
  mutate(time_value = time_value - tvmod, tvmod = NULL)
data_finalized <- inner_join(
  data_finalized |> select(-deaths),
  nchs_daily,
  by = c("geo_value", "time_value")
) |>
  arrange(time_value)

# Fix negatives and weekly reporting
versions <- sort(unique(data_archive$DT$version))
data_archive_fixed <- epix_slide(
  data_archive,
  .f = function(x, g, v) {
    mutate(
      x,
      deaths = fix_weekly_and_negatives(replace_na(deaths, 0), 1, 28)
    )
  },
  .versions = versions
) |>
  as_epi_archive()


makeExperimentRegistry(
  file.dir = "real-data-reg",
  packages = c("tidyverse", "epiprocess"),
  source = c(
    here::here("Code", "estimate.R"),
    here::here("Code", "helper_functions.R"),
    here::here("Code", "funs.R"),
    here::here("Code", "batch-funs.R")
  )
)


data_list <- enlist(
  raw_data = data_finalized |> select(-deaths),
  raw_deaths = data_finalized |> select(-hosps),
  data_archive = data_archive_fixed$DT |> select(-deaths) |> as_epi_archive(),
  archive_deaths = data_archive_fixed$DT |> select(-hosps) |> as_epi_archive(),
  meta_data,
  delay_mean_modifier = 0,
  delay_sds = 0.9,
  ws,
  d
)

create_real_problem_retro <- function(job, data, ...) {
  data$raw_deaths |>
    rename(sim_deaths = deaths) |>
    arrange(time_value) |>
    slice(-(1:data$d))
}

create_real_problem_rt <- function(job, data, ...) {
  data$archive_deaths$DT |>
    as.data.frame() |>
    mutate(sim_deaths = deaths) |> # algo removes `deaths`, then renames
    arrange(version, time_value) |>
    as_epi_archive()
}

addProblem("retro", data = data_list, fun = create_real_problem_retro)
addProblem("realtime", data = data_list, fun = create_real_problem_rt)
addAlgorithm("tf_retro", fun = trendfilter_retro_algo)
addAlgorithm("bench_retro", fun = benchmark_retro_algos)
addAlgorithm("tf_rt", fun = trendfilter_rt_algo)
addAlgorithm("bench_rt", fun = benchmark_rt_algos)
addExperiments(
  prob.designs = list(retro = tibble(region = region)),
  algo.designs = list(
    bench_retro = tibble(mod_idx = 1),
    tf_retro = tibble(mod_idx = 1, order = 0:2)
  )
)
addExperiments(
  prob.designs = list(realtime = tibble(region = region)),
  algo.designs = list(
    bench_rt = tibble(mod_idx = 1),
    tf_rt = tibble(mod_idx = 1, order = 0:2)
  )
)
ids <- getJobPars() |> unwrap()
ids$chunk <- rep(c(1, 1:3), 2)
# submitJobs(ids, resources = list(walltime = "12:0:0", memory = "4gb", ncpus = 1))
