library(tidyverse)
# replace_inf(), used so the MAE matches the main-text collection pipelines.
source(here::here("Code", "Analysis", "collection-funs.R"))

# MAE evaluation windows, matched to the stationary main-text pipelines:
#   retro: burn 60 days off each end (collect-retro.R), then restrict to dates
#          common to all methods (order-2 trend filtering fans out at the
#          unburned boundaries, which otherwise dominates its MAE).
#   rt:    drop the first 3 estimates only (grab_rt_maes, trim_left_tail = 3L);
#          the right edge is "now" and must be kept.
retro_burn <- 60L
rt_trim_left <- 3L

# ---- compute MAEs from time-varying experiment curves ----
res <- read_rds(here::here("Data", "batch-results", "time-var-delay.rds"))
gt <- read_rds(here::here("Data", "raw-data", "daily-hfrs.rds")) |>
  filter(geo_value %in% c("ca", "la", "pa", "ut")) |>
  rename(truth = hfr)

res <- res |>
  mutate(
    method = case_when(
      algorithm %in% c("tf-rt", "tf-retro") ~ paste0("tf-o", order),
      estimator == "conv" ~ "conv",
      estimator == "lagged" ~ "lagged",
      TRUE ~ paste(algorithm, estimator, sep = "-")
    )
  )

tvar_long <- res |>
  select(problem, region, method, spec, repl, res) |>
  unnest(res) |>
  rename(geo_value = region)

# Retrospective: burn 60 days off each end of every curve, then keep only the
# dates present for all methods (each curve = one geo x spec x repl x method).
retro_curves <- tvar_long |>
  filter(problem == "retro") |>
  arrange(geo_value, spec, repl, method, time_value) |>
  slice(-(1:retro_burn), -((n() - retro_burn):n()),
        .by = c(geo_value, spec, repl, method)) |>
  group_by(geo_value, spec, repl) |>
  mutate(.n_methods = n_distinct(method)) |>
  group_by(geo_value, spec, repl, time_value) |>
  filter(n() == .n_methods[1]) |>
  ungroup() |>
  select(-.n_methods)

# Real-time: drop only the first 3 estimates per curve (dates already aligned).
rt_curves <- tvar_long |>
  filter(problem == "rt") |>
  arrange(geo_value, spec, repl, method, time_value) |>
  slice(-(1:rt_trim_left), .by = c(geo_value, spec, repl, method))

tvar_curves <- bind_rows(retro_curves, rt_curves) |>
  inner_join(gt, by = c("geo_value", "time_value")) |>
  mutate(ae = abs(replace_inf(est) - truth))

tvar_maes <- tvar_curves |>
  summarise(mae = mean(ae, na.rm = TRUE),
            .by = c(problem, geo_value, method, spec, repl))

# ---- pull stationary MAEs for the same 4 states ----
load_stat <- function(file, prob, spec) {
  read_rds(here::here("Data", "batch-results", file)) |>
    as_tibble() |>
    filter(geo_value %in% c("ca", "la", "pa", "ut"),
           estimator != "tf_tuned_order") |>
    mutate(
      method = case_when(
        estimator == "trendfilter" ~ paste0("tf-o", order),
        TRUE ~ estimator
      ),
      problem = prob,
      spec = spec
    ) |>
    select(problem, geo_value, method, spec, mae)
}
stat_maes <- bind_rows(
  load_stat("retro-maes-mis.rds", "retro", "mis"),
  load_stat("retro-maes-ws.rds", "retro", "ws"),
  load_stat("rt-maes-mis.rds", "rt", "mis"),
  load_stat("rt-maes-ws.rds", "rt", "ws")
)

# ---- join and summarize ----
summary_tab <- bind_rows(
  tvar_maes |> mutate(experiment = "tvar"),
  stat_maes |> mutate(experiment = "stat")
) |>
  summarise(mean_mae = mean(mae),
            sd_mae = sd(mae),
            n = n(),
            .by = c(experiment, problem, geo_value, method, spec))

method_order <- c("tf-o0", "tf-o1", "tf-o2", "conv", "lagged")

# Wide print: rows = state x method, columns = experiment x spec
wide <- summary_tab |>
  mutate(method = factor(method, levels = method_order)) |>
  mutate(cell = sprintf("%.3f", mean_mae)) |>
  select(problem, geo_value, method, spec, experiment, cell) |>
  pivot_wider(names_from = c(experiment, spec),
              values_from = cell,
              names_glue = "{experiment}_{spec}") |>
  arrange(problem, geo_value, method) |>
  select(problem, geo_value, method,
         stat_mis, tvar_mis, stat_ws, tvar_ws)

cat("\n=== Mean MAE by state x method (rounded to 3 d.p.) ===\n")
cat("    stat_*: stationary-delay experiment (existing)\n")
cat("    tvar_*: time-varying-delay experiment (new)\n\n")
for (prob in c("retro", "rt")) {
  cat("---- problem:", prob, "----\n")
  print(wide |> filter(problem == prob) |> select(-problem), n = Inf)
  cat("\n")
}

# ---- compact per-state summary: avg MAE across methods ----
cat("\n=== Per-state averages across methods (mean MAE) ===\n")
print(
  summary_tab |>
    summarise(mean_mae = mean(mean_mae),
              .by = c(experiment, problem, geo_value, spec)) |>
    mutate(cell = sprintf("%.3f", mean_mae)) |>
    select(-mean_mae) |>
    pivot_wider(names_from = c(experiment, spec), values_from = cell,
                names_glue = "{experiment}_{spec}") |>
    arrange(problem, geo_value),
  n = Inf
)

# ---- save the tidy tables ----
out_dir <- here::here("Data", "batch-results")
write_rds(tvar_maes, file.path(out_dir, "time-var-maes.rds"))
write_rds(summary_tab, file.path(out_dir, "time-var-mae-summary.rds"))
cat("\nwrote time-var-maes.rds and time-var-mae-summary.rds\n")
