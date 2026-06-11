library(batchtools)
library(tidyverse)
library(epiprocess)
library(data.table)
source(here::here("Code", "helper_functions.R"))
source(here::here("Code", "estimate.R"))
source(here::here("Code", "funs.R"))
source(here::here("Code", "batch-funs.R"))

regions <- c("ca", "pa", "la", "ut")
variants <- c("Original", "Alpha", "Delta", "Omicron")
d <- 60L
ws <- c(1, 7, 15, 21, 29)
ord <- 0:2


# produced by the commented code below
ca_voc_meta <- tibble(
  variant = variants,
  retro_lag_argmax = c(8, 5, 13, 9),
  rt_lag_1se = c(18, 6, 22, 21)
) |>
  rowwise() |>
  mutate(
    retro_delay = list(make_delay_distr(
      retro_lag_argmax,
      0.9 * retro_lag_argmax,
      d
    )),
    rt_delay = list(make_delay_distr(rt_lag_1se, 0.9 * rt_lag_1se, d))
  ) |>
  ungroup()

# finalized <- read_rds(here::here("Data", "raw-data", "finalized.rds"))
# rt_cleaned <- read_rds(here::here("Data", "raw-data", "cleaned-latest.rds")) |>
#   ungroup()

# variant_props <- read_rds(here::here(
#   "Data",
#   "raw-data",
#   "daily-variant-props.rds"
# )) |>
#   pivot_longer(Original:Omicron, names_to = "variant", values_to = "prop") |>
#   summarise(dom_voc = variant[which.max(prop)], .by = c(time_value, geo_value))

# retro_meta_data <- finalized |>
#   left_join(
#     variant_props,
#     by = join_by(geo_value, time_value)
#   ) |>
#   summarise(
#     retro_lag_argmax = pmin(pmax(5, compute_max_acf_lag(hosps, deaths)), 35),
#     .by = c(geo_value, dom_voc)
#   ) |>
#   filter(geo_value %in% regions)

# rt_meta_data <- rt_cleaned |>
#   left_join(
#     variant_props,
#     by = join_by(geo_value, time_value)
#   ) |>
#   summarise(
#     dispersion = summary(
#       glm(
#         deaths_fix_reporting ~ deaths_mu + as.factor(wday(time_value)),
#         family = quasipoisson()
#       )
#     )$dispersion,
#     rt_lag_argmax = pmin(pmax(6, compute_max_acf_lag(hosps, deaths_7dav)), 35),
#     rt_lag_1se = pmin(pmax(6, compute_1se_acf_lag(hosps, deaths_7dav)), 35),
#     .by = c(geo_value, dom_voc)
#   ) |>
#   filter(geo_value %in% regions)

meta_data <- read_rds(here::here("Data", "raw-data", "all-meta-data.rds"))

# Create time varying delays
variant_props <- read_rds(here::here(
  "Data",
  "raw-data",
  "daily-variant-props.rds"
)) |>
  pivot_longer(Original:Omicron, names_to = "variant", values_to = "prop") |>
  filter(geo_value %in% regions) |>
  left_join(
    ca_voc_meta |> select(variant, retro_delay, rt_delay),
    by = "variant"
  ) |>
  summarise(
    retro_delay = list(colSums(prop * do.call(rbind, retro_delay))),
    rt_delay = list(colSums(prop * do.call(rbind, rt_delay))),
    .by = c(geo_value, time_value)
  )

gt_hfrs <- read_rds(here::here("Data", "raw-data", "daily-hfrs.rds")) |>
  filter(geo_value %in% regions)
raw_data <- read_rds(here::here("Data", "raw-data", "finalized.rds")) |>
  filter(geo_value %in% regions) |>
  select(-deaths)
data_archive <- read_rds(here::here("Data", "raw-data", "full-archive.rds")) |>
  filter(geo_value %in% regions) |>
  as_epi_archive()
data_archive$geo_type <- "nation"


data_list <- enlist(
  raw_data,
  data_archive,
  gt_hfrs,
  meta_data,
  variant_props,
  ws,
  d
)

pd_retro <- list(retro = data.frame(region = regions, time_var = TRUE))
pd_rt <- list(rt = data.frame(region = regions, time_var = TRUE))
ad_retro <- list(
  `tf-retro` = data.table::CJ(order = c(0, 1, 2), spec = c("mis", "ws")),
  `bench-retro` = data.frame(spec = c("mis", "ws"))
)

ad_rt <- list(
  `tf-rt` = data.table::CJ(order = c(0, 1, 2), spec = c("mis", "ws")),
  `bench-rt` = data.frame(spec = c("mis", "ws"))
)


# Set up the experimental design
makeExperimentRegistry(
  "severity-estimation-tvar-delay",
  source = c(
    here::here("Code", "helper_functions.R"),
    here::here("Code", "estimate.R"),
    here::here("Code", "funs.R"),
    here::here("Code", "batch-funs.R")
  ),
  packages = c("tidyverse", "epiprocess")
)
addProblem(
  "retro",
  fun = create_regional_problem_retro,
  data = data_list,
  seed = 123
)
addProblem(
  "rt",
  fun = create_regional_problem_rt,
  data = data_list,
  seed = 123
)

addAlgorithm("tf-rt", fun = trendfilter_rt_algo)
addAlgorithm("bench-rt", fun = benchmark_rt_algos)
addAlgorithm("tf-retro", fun = trendfilter_retro_algo)
addAlgorithm("bench-retro", fun = benchmark_retro_algos)
addExperiments(pd_rt, ad_rt, repls = 10L)
addExperiments(pd_retro, ad_retro, repls = 10L)
