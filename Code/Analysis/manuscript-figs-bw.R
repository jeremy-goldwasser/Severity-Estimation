# Setup -------------------------------------------------------------------
# Black-and-white version of manuscript-figs.R. Series are distinguished by
# greyscale shade combined with linetype (lines) or point shape (scatter)
# instead of colour. Output paths are identical to the colour script, so
# running this populates B&W figures where the colour ones used to live (the
# colour originals were moved to Figures/Color/).
library(tidyverse)
library(cowplot)
library(epidatr)
library(usmap)
library(glue)
library(purrr)
library(readr)
library(RcppRoll)
library(patchwork)
data("statepop")
# Greyscale shades. Co-occurring series get distinct shades; where two share a
# shade (Ground Truth vs Deconv-0) the linetype below keeps them apart.
colours <- c(
  `Ground Truth` = "grey0",
  `Conv Ratio` = "grey55",
  `Lagged Ratio` = "grey75",
  `Deconv-0` = "grey0",
  `Deconv-1` = "grey35",
  `Deconv-2` = "grey60"
)
linetypes <- c(
  `Ground Truth` = "solid",
  `Conv Ratio` = "dashed",
  `Lagged Ratio` = "dotted",
  `Deconv-0` = "dotdash",
  `Deconv-1` = "longdash",
  `Deconv-2` = "twodash"
)
shapes <- c(
  `Ground Truth` = 15,
  `Conv Ratio` = 16,
  `Lagged Ratio` = 17,
  `Deconv-0` = 18,
  `Deconv-1` = 8,
  `Deconv-2` = 4
)
methods <- names(colours)[-1]

# Tables print LaTeX via knitr::kable when available; this script is about
# figures, so fall back to a no-op if knitr is not installed.
if (requireNamespace("knitr", quietly = TRUE)) {
  kable <- knitr::kable
} else {
  kable <- function(...) invisible(NULL)
}

name_map <- c(
  "Deconv-0" = "0",
  "Deconv-1" = "1",
  "Deconv-2" = "2",
  "Deconv-Tuned" = "Tuned",
  "Conv Ratio" = "Conv Ratio",
  "Lagged Ratio" = "Lagged Ratio"
)


# Potentially required package installs. Not all may be used
# [1] "batchtools"   "cowplot"      "CVXR"         "data.table"   "dplyr"        "epidatr"
# [7] "epiprocess"   "extraDistr"   "forecast"     "future.apply" "genlasso"     "ggplot2"
# [13] "grid"         "gridExtra"    "here"         "lubridate"    "Matrix"       "parallel"
# [19] "patchwork"    "purrr"        "rmarkdown"    "scales"       "splines"      "stats"
# [25] "tibble"       "tidyr"        "tidyverse"    "usmap"        "VGAM"         "zoo"

source(here::here("Code", "helper_functions.R"))
source(here::here("Code", "funs.R"))
source(here::here("Code", "batch-funs.R"))
source(here::here("Code", "estimate.R"))
d <- 60


retro_data <- read_rds(here::here("Data", "raw-data", "finalized.rds"))
cleaned_data <- read_rds(here::here("Data", "raw-data", "cleaned-latest.rds"))
rt_data <- read_rds(here::here("Data", "raw-data", "full-archive.rds")) |>
  as_epi_archive()
gt_hfrs <- read_rds(here::here("Data", "raw-data", "daily-hfrs.rds"))
meta_data <- read_rds(here::here("Data", "raw-data", "all-meta-data.rds"))

retro_curves_og <- read_rds(here::here(
  "Data",
  "batch-results",
  "retro-example-curves.rds"
))


# Figures: Real v simulated deaths --------------------------------------

regions <- c("ca", "pa", "la", "ut")
region_titles <- c("California", "Pennsylvania", "Louisiana", "Utah")

retro_sim <- map(regions, \(r) {
  create_regional_problem_retro(
    enlist(raw_data = retro_data, gt_hfrs, meta_data, d = 60), # added d
    r
  )
}) |>
  bind_rows()

rt_sim <- map(regions, \(r) {
  create_regional_problem_rt(
    enlist(data_archive = rt_data, gt_hfrs, meta_data, d = 60),
    r
  ) |>
    epix_as_of_current()
}) |>
  bind_rows()

names(region_titles) <- regions

# Add noise to the constant, weekly data
nchs_deaths_daily <- retro_data |>
  filter(wday(time_value) == 7, !is.na(deaths)) |>
  rowwise() |>
  mutate(
    daily = list(rmultinom(1, deaths * 7, prob = rep(1, 7))[, 1]),
    days = list(seq(time_value - 6, time_value, by = "1 day"))
  ) |>
  unnest(c(daily, days)) %>%
  transmute(
    geo_value,
    time_value = days,
    deaths = daily
  ) |>
  ungroup()

inner_join(nchs_deaths_daily, retro_sim, by = join_by(time_value, geo_value)) |>
  pivot_longer(c(deaths, sim_deaths)) |>
  mutate(
    geo_value = factor(geo_value, levels = regions),
    name = case_when(name == "deaths" ~ "Real", TRUE ~ "Simulated")
  ) |>
  ggplot(aes(time_value, value, color = name, shape = name)) +
  geom_point(alpha = .6) +
  facet_wrap(
    ~geo_value,
    nrow = 2,
    scales = "free_y",
    labeller = labeller(geo_value = region_titles)
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.05))) +
  theme_bw() +
  xlab("") +
  ylab("Deaths") +
  theme(strip.background = element_blank()) +
  scale_color_manual(values = c(Real = "grey20", Simulated = "grey65"), name = "") +
  scale_shape_manual(values = c(Real = 16, Simulated = 17), name = "") +
  guides(
    colour = guide_legend(override.aes = list(alpha = 1, size = 3)),
    shape = guide_legend(override.aes = list(alpha = 1, size = 3))
  )

ggsave(
  here::here("Figures", "BW", "Retrospective", "sim_vs_real_deaths.pdf"),
  width = 8,
  height = 5
)

cleaned_data |>
  select(time_value, geo_value, deaths = deaths_fix_reporting) |>
  inner_join(rt_sim, by = join_by(time_value, geo_value)) |>
  pivot_longer(c(deaths, sim_deaths)) |>
  mutate(
    geo_value = factor(geo_value, levels = regions),
    name = case_when(name == "deaths" ~ "Real", TRUE ~ "Simulated")
  ) |>
  ggplot(aes(time_value, value, color = name, shape = name)) +
  geom_point(alpha = .6) +
  facet_wrap(
    ~geo_value,
    nrow = 2,
    scales = "free_y",
    labeller = labeller(geo_value = region_titles)
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.05))) +
  theme_bw() +
  xlab("") +
  ylab("Deaths") +
  theme(strip.background = element_blank()) +
  scale_color_manual(values = c(Real = "grey20", Simulated = "grey65"), name = "") +
  scale_shape_manual(values = c(Real = 16, Simulated = 17), name = "") +
  guides(
    colour = guide_legend(override.aes = list(alpha = 1, size = 3)),
    shape = guide_legend(override.aes = list(alpha = 1, size = 3))
  )

# Note the warning: some negative observed deaths are below the 0-line

ggsave(
  here::here("Figures", "BW", "Real_Time", "sim_vs_real_deaths.pdf"),
  width = 8,
  height = 5
)

# Table 1 & Figure - Retrospective severity estimator comparison ------------------
# ---- Load Data ----
retro_maes <- read_rds(here::here(
  "Data",
  "batch-results",
  "retro-maes-ws.rds"
)) |>
  mutate(
    est_name = case_when(
      is.na(order) ~ estimator,
      order == 0 ~ "Deconv-0",
      order == 1 ~ "Deconv-1",
      order == 2 ~ "Deconv-2"
    )
  ) |>
  mutate(
    est_name = case_match(
      est_name,
      "conv" ~ "Conv Ratio",
      "lagged" ~ "Lagged Ratio",
      .default = est_name
    )
  )

# add tuned order
retro_maes <- retro_maes |>
  bind_rows(
    retro_maes |>
      filter(estimator == "trendfilter") |>
      group_by(geo_value, repl) |>
      summarize(
        mae = mae[cverr == min(cverr)],
        est_name = "Deconv-Tuned",
        estimator = "trendfilter",
        .groups = "drop"
      )
  )

# Oracle version
retro_maes_oracle <- read_rds(here::here(
  "Data",
  "batch-results",
  "retro-maes-ws-oracle.rds"
)) |>
  mutate(
    est_name = case_when(
      is.na(order) ~ estimator,
      order == 0 ~ "Deconv-0",
      order == 1 ~ "Deconv-1",
      order == 2 ~ "Deconv-2"
    )
  ) |>
  mutate(
    est_name = case_match(
      est_name,
      "conv" ~ "Conv Ratio",
      "lagged" ~ "Lagged Ratio",
      .default = est_name
    )
  )

retro_maes_oracle <- retro_maes_oracle |>
  bind_rows(
    retro_maes_oracle |>
      filter(estimator == "trendfilter") |>
      group_by(geo_value, repl) |>
      summarize(
        mae = mae[cverr == min(cverr)],
        est_name = "Deconv-Tuned",
        estimator = "trendfilter",
        .groups = "drop"
      )
  ) |>
  mutate(Tuning = "Oracle")

# ---- Combine ----
retro_maes_all <- bind_rows(retro_maes, retro_maes_oracle) |>
  mutate(Tuning = if_else(is.na(Tuning), "CV", Tuning))

# ---- Estimators of interest ----
all_estimators <- c(
  "Lagged Ratio",
  "Conv Ratio",
  "Deconv-0",
  "Deconv-1",
  "Deconv-2",
  "Deconv-Tuned"
)

# ---- Function to make LaTeX table ----
make_table_data <- function(df) {
  # MAE row
  tab_maes <- df |>
    group_by(est_name, geo_value) |>
    summarize(
      mean_mae_region = mean(mae, na.rm = TRUE),
      var_mae_region = var(mae, na.rm = TRUE),
      n_region = sum(!is.na(mae)),
      .groups = "drop"
    ) |>
    group_by(est_name) |>
    summarize(
      mean_mae = mean(mean_mae_region, na.rm = TRUE),
      Vxbar = sum(n_region * var_mae_region, na.rm = TRUE) / (sum(n_region)^2),
      SE = sqrt(Vxbar),
      .groups = "drop"
    ) |>
    mutate(across(where(is.numeric), ~ .x * 1000)) |>
    mutate(cell = sprintf("%.1f ± %.1f", mean_mae, SE)) |>
    select(est_name, cell) |>
    pivot_wider(names_from = est_name, values_from = cell) |>
    mutate(row = "MAE × 10^3") |>
    relocate(row)

  # Improvements
  tab_imps <- purrr::map_dfr(
    c("Conv Ratio", "Lagged Ratio"),
    function(bench_name) {
      df |>
        group_by(geo_value, repl) |>
        mutate(
          mae_bench = mae[est_name == bench_name],
          I = (mae_bench - mae) / mae_bench * 100
        ) |>
        ungroup() |>
        filter(est_name %in% all_estimators) |>
        group_by(est_name, geo_value) |>
        summarize(
          var_I_region = var(I, na.rm = TRUE),
          mean_I_region = mean(I, na.rm = TRUE),
          n_region = sum(!is.na(I)),
          .groups = "drop"
        ) |>
        group_by(est_name) |>
        summarize(
          mean_I = mean(mean_I_region, na.rm = TRUE),
          Vxbar = sum(n_region * var_I_region, na.rm = TRUE) /
            (sum(n_region)^2),
          SE = sqrt(Vxbar),
          .groups = "drop"
        ) |>
        mutate(
          cell = sprintf("%.1f ± %.1f", mean_I, SE),
          row = if_else(
            bench_name == "Conv Ratio",
            "Improvement over CR (%)",
            "Improvement over LR (%)"
          )
        ) |>
        select(est_name, cell, row)
    }
  ) |>
    pivot_wider(names_from = est_name, values_from = cell) |>
    relocate(row)

  # Combine tibble
  tab <- bind_rows(tab_maes, tab_imps) |>
    select(row, all_of(all_estimators))

  tab
}

# --- Then format with kable only when needed ---
make_table <- function(df, tuning_label) {
  tab <- make_table_data(df)
  kable(
    tab,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    caption = paste("Retrospective estimators,", tuning_label, "tuning")
  )
}

tab_cv <- retro_maes_all %>% filter(Tuning == "CV")
tab_oracle <- retro_maes_all %>% filter(Tuning == "Oracle")

# Inspect in R as a tibble
make_table_data(tab_cv) # Table 1
make_table_data(tab_oracle) # Table 3

# Get LaTeX when ready
make_table(tab_cv, "CV")
make_table(tab_oracle, "Oracle")

# Figure
pop <- statepop |>
  select(geo_value = abbr, pop = pop_2022) |>
  mutate(geo_value = tolower(geo_value)) |>
  add_row(geo_value = "us", pop = 33 * 1e7)

retro_maes_by_pop <- retro_maes |>
  select(mae, geo_value, est_name) |>
  summarise(mae = mean(mae, na.rm = TRUE), .by = c(geo_value, est_name)) |>
  left_join(pop, by = "geo_value")

ylim <- range(retro_maes_by_pop$mae)
retro_maes_by_pop <- retro_maes_by_pop |>
  mutate(est_name = fct_relevel(est_name, methods))

pPop1 <- ggplot(
  retro_maes_by_pop |>
    filter(est_name %in% c("Conv Ratio", "Lagged Ratio", "Deconv-0")),
  aes(x = pop, y = mae, color = est_name, shape = est_name)
) +
  # geom_smooth(se = FALSE, method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(
    breaks = c("Conv Ratio", "Lagged Ratio", "Deconv-0"),
    values = colours
  ) +
  scale_shape_manual(
    breaks = c("Conv Ratio", "Lagged Ratio", "Deconv-0"),
    values = shapes
  ) +
  labs(x = "Population", y = "MAE", color = "", shape = "") +
  theme_bw() +
  scale_y_continuous(limits = c(0, ylim[2]), expand = expansion(c(0, 0.05))) +
  scale_x_log10(labels = scales::label_log()) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

pPop2 <- ggplot(
  retro_maes_by_pop |>
    filter(est_name %in% c("Deconv-0", "Deconv-1", "Deconv-2")), #, "Deconv-Tuned"
  aes(x = pop, y = mae, color = est_name, shape = est_name)
) +
  # geom_smooth(se = FALSE, method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(values = colours) +
  scale_shape_manual(values = shapes) +
  labs(x = "Population", y = "MAE", color = "", shape = "") +
  theme_bw() +
  scale_y_continuous(limits = c(0, ylim[2]), expand = expansion(c(0, 0.05))) +
  scale_x_log10(labels = scales::label_log()) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

plot_grid(
  pPop1,
  pPop2,
  nrow = 1,
  rel_widths = c(0.5, 0.5),
  align = "h",
  labels = "AUTO"
)
ggsave(
  here::here("Figures", "BW", "Retrospective", "improvements.pdf"),
  width = 8,
  height = 4
)

# Figure - Retrospective rates  ------------------------------------

single_region <- "pa"
region_curves1 <- filter(retro_curves_og, region == single_region) |>
  mutate(
    est_name = case_when(
      estimator == "conv" ~ "Conv Ratio",
      estimator == "lagged" ~ "Lagged Ratio",
      order == 0 ~ "Deconv-0",
      order == 1 ~ "Deconv-1",
      order == 2 ~ "Deconv-2"
    ),
    rl = case_when(
      estimator %in% c("conv", "lagged") ~ "l",
      algorithm == "trendfilter" & order == 0 ~ "b",
      algorithm == "trendfilter" ~ "r",
      TRUE ~ "l"
    )
  ) |>
  # Select single replicate
  filter(repl==1)

region_curves2 <- region_curves1 |>
  unnest(result) |>
  select(time_value, est, est_name, rl) |>
  bind_rows(
    gt_hfrs |>
      filter(geo_value == single_region) |>
      mutate(est_name = "Ground Truth", rl = "l") |>
      rename(est = hfr) |>
      select(-geo_value)
  ) |>
  bind_rows(
    gt_hfrs |>
      filter(geo_value == single_region) |>
      mutate(est_name = "Ground Truth", rl = "r") |>
      rename(est = hfr) |>
      select(-geo_value)
  ) |>
  mutate(
    est_name = fct_relevel(
      as.factor(est_name),
      "Ground Truth",
      "Conv Ratio",
      "Lagged Ratio",
      # "Deconv-Tuned",
      "Deconv-0",
      "Deconv-1",
      "Deconv-2"
    ),
    ww = ifelse(est_name == "Ground Truth", "a", "b")
  )

region_curves3 <- region_curves2 |>
  bind_rows(region_curves2 |> filter(rl == "b") |> mutate(rl = "l")) |>
  bind_rows(region_curves2 |> filter(rl == "b") |> mutate(rl = "r")) |>
  filter(rl != "b")

bounds <- region_curves3 |>
  filter(est_name == "Conv Ratio") |>
 summarize(xmin = min(time_value), xmax = max(time_value))

region_curves <- region_curves3 |> drop_na() |>
  filter(time_value >= bounds$xmin, time_value <= bounds$xmax)
ylim <- range(region_curves3$est)

region_curves |> group_by(est_name) |>
  summarize(min(time_value), max(time_value))

p1 <- ggplot(
  region_curves |> filter(rl == "l") %>%
    arrange(time_value) %>%
    filter(as.integer(time_value - min(time_value)) %% 7 == 0),
  aes(time_value, est, color = est_name, alpha = est_name, linetype = est_name)
) +
  geom_line(linewidth = 0.7) +
  scale_alpha_manual(
    breaks = names(colours),
    values = c(1, rep(0.8, 10)),
    guide = "none"
  ) +
  scale_color_manual(
    breaks = c("Conv Ratio", "Lagged Ratio", "Deconv-0"),
    values = colours
  ) +
  scale_linetype_manual(
    breaks = c("Conv Ratio", "Lagged Ratio", "Deconv-0"),
    values = linetypes
  ) +
  theme_bw() +
  ylab("HFR") +
  xlab("") +
  scale_y_continuous(limits = ylim) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(t = -10),
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

p2 <- ggplot(
  region_curves |> filter(rl == "r", est_name != "Deconv-Tuned"),
  aes(time_value, est, color = est_name, alpha = est_name, linetype = est_name)
) +
  geom_line(linewidth = 0.7) +
  scale_alpha_manual(
    breaks = names(colours),
    values = c(1, rep(0.8, 10)),
    guide = "none"
  ) +
  scale_color_manual(
    breaks = c("Deconv-0", "Deconv-1", "Deconv-2"),
    values = colours
  ) +
  scale_linetype_manual(
    breaks = c("Deconv-0", "Deconv-1", "Deconv-2"),
    values = linetypes
  ) +
  theme_bw() +
  ylab("HFR") +
  xlab("") +
  scale_y_continuous(limits = ylim) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(t = -10)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

plot_grid(
  p1,
  p2,
  nrow = 1,
  rel_widths = c(0.5, 0.5),
  align = "h",
  labels = "AUTO"
)
ggsave(
  here::here("Figures", "BW", "Retrospective", "single_curve.pdf"),
  width = 8,
  height = 4
)

# tf_dates <- region_curves |> filter(est_name=="Deconv-0") |> pull(time_value)
cr_dates <- region_curves |>
  filter(est_name == "Conv Ratio") |>
  pull(time_value) |>
  unique()
region_gt <- region_curves |>
  filter(est_name == "Ground Truth", time_value %in% cr_dates, rl == "l") |>
  pull(est)
region_tf0 <- region_curves |>
  filter(est_name == "Deconv-0", time_value %in% cr_dates) |>
  pull(est)
region_tf1 <- region_curves |>
  filter(est_name == "Deconv-1", time_value %in% cr_dates) |>
  pull(est)
region_tf2 <- region_curves |>
  filter(est_name == "Deconv-2", time_value %in% cr_dates) |>
  pull(est)
region_tft <- region_curves |>
  filter(est_name == "Deconv-Tuned", time_value %in% cr_dates) |>
  pull(est)
region_cr <- region_curves |>
  filter(est_name == "Conv Ratio", time_value %in% cr_dates) |>
  pull(est)
mae(region_tf0, region_gt)
mae(region_tf1, region_gt)
mae(region_tf2, region_gt)
mae(region_tft, region_gt)
mae(region_cr, region_gt)

# Figure - Real time severity estimator comparison --------------------------

rt_maes <- read_rds(here::here("Data", "batch-results", "rt-maes-ws.rds")) |>
  mutate(
    est_name = case_when(
      order == 0 ~ "Deconv-0",
      order == 1 ~ "Deconv-1",
      order == 2 ~ "Deconv-2",
      estimator == "tf_tuned_order" ~ "Deconv-Tuned",
      is.na(order) ~ estimator,
    )
  )
conv_mae <- rt_maes |>
  filter(estimator == "conv") |>
  select(conv_mae = mae, geo_value, repl)
lagged_mae <- rt_maes |>
  filter(estimator == "lagged") |>
  select(lagged_mae = mae, geo_value, repl)
rt_plot_df <- bind_rows(rt_maes) |>
  filter(estimator %in% c("trendfilter", "tf_tuned_order")) |>
  left_join(conv_mae, by = join_by(geo_value, repl)) |>
  left_join(lagged_mae, by = join_by(geo_value, repl))

# Table 2: Real-time MAEs and percent improvements
rt_long <- rt_plot_df |>
  tidyr::pivot_longer(
    cols = c(conv_mae, lagged_mae),
    names_to = "baseline",
    values_to = "mae_val"
  ) |>
  mutate(
    est_name = dplyr::case_match(
      baseline,
      "conv_mae" ~ "Conv Ratio",
      "lagged_mae" ~ "Lagged Ratio",
      .default = baseline
    )
  ) |>
  select(job.id, mae = mae_val, geo_value, order, repl, estimator, est_name)

rt_all <- bind_rows(
  rt_plot_df |>
    select(job.id, mae, geo_value, order, repl, estimator, est_name),
  rt_long
)

all_estimators <- c(
  "Lagged Ratio",
  "Conv Ratio",
  "Deconv-0",
  "Deconv-1",
  "Deconv-2",
  "Deconv-Tuned"
)

# --- MAE row ---
tab_maes <- rt_all |>
  filter(est_name %in% all_estimators) |>
  group_by(est_name, geo_value) |>
  summarize(
    mean_mae_region = mean(mae, na.rm = TRUE),
    var_mae_region = var(mae, na.rm = TRUE),
    n_region = sum(!is.na(mae)),
    .groups = "drop"
  ) |>
  group_by(est_name) |>
  summarize(
    mean_mae = mean(mean_mae_region, na.rm = TRUE),
    Vxbar = sum(n_region * var_mae_region, na.rm = TRUE) / (sum(n_region)^2),
    SE = sqrt(Vxbar),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), ~ .x * 1000)) |>
  mutate(cell = sprintf("%.1f $\\pm$ %.1f", mean_mae, SE)) |>
  select(est_name, cell) |>
  tidyr::pivot_wider(names_from = est_name, values_from = cell) |>
  mutate(row = "MAE $\\times 10^3$") |>
  relocate(row)

# --- Improvement rows (TF + baselines) ---
benchmarks <- c("Conv Ratio", "Lagged Ratio")

tab_imps <- purrr::map_dfr(benchmarks, function(bench_name) {
  bench_df <- rt_all |>
    filter(est_name == bench_name) |>
    distinct(geo_value, repl, .keep_all = TRUE) |>
    transmute(geo_value, repl, mae_bench = mae, benchmark = est_name)

  rt_all |>
    filter(est_name %in% all_estimators) |>
    distinct(geo_value, repl, est_name, mae, .keep_all = TRUE) |>
    mutate(benchmark = bench_name) |>
    left_join(bench_df, by = c("geo_value", "repl", "benchmark")) |>
    mutate(I = (mae_bench - mae) / mae_bench * 100) |>
    group_by(est_name, geo_value, benchmark) |>
    summarize(
      var_I_region = var(I, na.rm = TRUE),
      mean_I_region = mean(I, na.rm = TRUE),
      n_region = sum(!is.na(I)),
      .groups = "drop"
    ) |>
    group_by(est_name, benchmark) |>
    summarize(
      mean_I = mean(mean_I_region, na.rm = TRUE),
      Vxbar = sum(n_region * var_I_region, na.rm = TRUE) / (sum(n_region)^2),
      SE = sqrt(Vxbar),
      .groups = "drop"
    ) |>
    mutate(
      cell = sprintf("%.1f $\\pm$ %.1f", mean_I, SE),
      row = if_else(
        bench_name == "Conv Ratio",
        "Improvement over CR (\\%)",
        "Improvement over LR (\\%)"
      )
    ) |>
    select(est_name, cell, row)
}) |>
  tidyr::pivot_wider(names_from = est_name, values_from = cell) |>
  relocate(row)

# --- Combine into final table ---
rt_table <- bind_rows(tab_maes, tab_imps) |>
  select(row, all_of(all_estimators))

# --- Output in LaTeX (and see in console first) ---
rt_table # Table 2
kable(
  rt_table,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Real-time estimators, CV tuning"
)

# Figure: Boxplots of improvements and MAE by population

rt_maes_by_pop <- rt_maes |>
  # filter(est_name %in% c("conv", "lagged", "Deconv-Tuned")) |>
  select(mae, geo_value, est_name) |>
  mutate(
    est_name = case_when(
      est_name == "conv" ~ "Conv Ratio",
      est_name == "lagged" ~ "Lagged Ratio",
      TRUE ~ est_name
      # est_name == "Deconv-Tuned" ~ "Deconv-Tuned"
    ),
    est_name = fct_relevel(as.factor(est_name), "Lagged Ratio")
  ) |>
  summarise(mae = median(mae, na.rm = TRUE), .by = c(geo_value, est_name)) |>
  left_join(pop, by = "geo_value")

ylim <- range(rt_maes_by_pop$mae)
rt_maes_by_pop <- rt_maes_by_pop |>
  mutate(est_name = str_replace(est_name, "^TF", "Deconv"))

pPop1 <- ggplot(
  rt_maes_by_pop |>
    filter(est_name %in% c("Conv Ratio", "Lagged Ratio", "Deconv-0")),
  aes(x = pop, y = mae, color = est_name, shape = est_name)
) +
  # geom_smooth(se = FALSE, method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(
    breaks = c("Conv Ratio", "Lagged Ratio", "Deconv-0"),
    values = colours
  ) +
  scale_shape_manual(
    breaks = c("Conv Ratio", "Lagged Ratio", "Deconv-0"),
    values = shapes
  ) +
  labs(x = "Population", y = "MAE", color = "", shape = "") +
  theme_bw() +
  scale_y_continuous(limits = c(0, ylim[2]), expand = expansion(c(0, 0.05))) +
  scale_x_log10(labels = scales::label_log()) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

pPop2 <- ggplot(
  rt_maes_by_pop |>
    filter(est_name %in% c("Deconv-0", "Deconv-1", "Deconv-2")),
  aes(x = pop, y = mae, color = est_name, shape = est_name)
) +
  # geom_smooth(se = FALSE, method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(values = colours) +
  scale_shape_manual(values = shapes) +
  labs(x = "Population", y = "MAE", color = "", shape = "") +
  theme_bw() +
  scale_y_continuous(limits = c(0, ylim[2]), expand = expansion(c(0, 0.05))) +
  scale_x_log10(labels = scales::label_log()) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

plot_grid(
  pPop1,
  pPop2,
  nrow = 1,
  rel_widths = c(0.5, 0.5),
  align = "h",
  labels = "AUTO"
)
ggsave(
  here::here("Figures", "BW", "Real_Time", "improvements.pdf"),
  width = 8,
  height = 4
)

# Figure - Real time rates  ----------------

# filter to only the nowcast date.
rt_curves <- read_rds(here::here(
  "Data",
  "batch-results",
  "rt-example-curves.rds"
))
rt_curves <- filter(rt_curves, repl == 1)
rt_curves <- filter(rt_curves, wday(time_value) == 1L)
single_region <- "pa" 
region_curves <- filter(rt_curves, region == single_region) |>
  mutate(
    est_name = case_when(
      estimator == "conv" ~ "Conv Ratio",
      estimator == "lagged" ~ "Lagged Ratio",
      estimator == "tf_tuned_order" ~ "Deconv-Tuned",
      order == 0 ~ "Deconv-0",
      order == 1 ~ "Deconv-1",
      order == 2 ~ "Deconv-2"
    ),
    # rl = case_when(estimator != "trendfilter" ~ "l", TRUE ~ "r")
    rl = case_when(
      estimator %in% c("conv", "lagged") ~ "l",
      #estimator == "tf_tuned_order" ~ "l",
      estimator == "trendfilter" & order == 0 ~ "b",
      estimator %in% c("tf_tuned_order", "trendfilter") ~ "r"
    )
  ) |>
  bind_rows(
    gt_hfrs |>
      filter(geo_value == single_region) |>
      mutate(est_name = "Ground Truth", rl = "l") |>
      rename(est = hfr)
  ) |>
  bind_rows(
    gt_hfrs |>
      filter(geo_value == single_region) |>
      mutate(est_name = "Ground Truth", rl = "r") |>
      rename(est = hfr)
  ) |>
  mutate(
    est_name = fct_relevel(
      as.factor(est_name),
      "Ground Truth",
      "Conv Ratio",
      "Lagged Ratio",
      "Deconv-Tuned",
      "Deconv-0",
      "Deconv-1",
      "Deconv-2"
    ),
    ww = ifelse(est_name == "Ground Truth", "a", "b")
  )
region_curves <- region_curves |>
  bind_rows(region_curves |> filter(rl == "b") |> mutate(rl = "l")) |>
  bind_rows(region_curves |> filter(rl == "b") |> mutate(rl = "r")) |>
  filter(rl != "b")

region_curves <- region_curves |>
  mutate(est_name = str_replace(est_name, "^TF", "Deconv")) |>
  select(time_value, est, est_name, rl, ww)

bounds <- region_curves |>
  filter(est_name == "Conv Ratio") |>
  summarize(xmin = min(time_value), xmax = max(time_value))

region_curves <- region_curves |> drop_na(est) |>
  filter(time_value >= bounds$xmin, time_value <= bounds$xmax)
ylim <- range(region_curves$est)

region_curves |> group_by(est_name) |>
  summarize(min(time_value), max(time_value))

p1 <- ggplot(
  region_curves |> filter(rl == "l"),
  aes(time_value, est, color = est_name, alpha = est_name, linetype = est_name)
) +
  geom_line(linewidth = 0.7) +
  scale_alpha_manual(
    breaks = names(colours),
    values = c(1, rep(0.8, 10)),
    guide = "none"
  ) +
  scale_color_manual(
    breaks = c("Conv Ratio", "Lagged Ratio", "Deconv-0"),
    values = colours
  ) +
  scale_linetype_manual(
    breaks = c("Conv Ratio", "Lagged Ratio", "Deconv-0"),
    values = linetypes
  ) +
  theme_bw() +
  ylab("HFR") +
  xlab("") +
  scale_y_continuous(limits = ylim) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.margin = margin(t = -10)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

p2 <- ggplot(
  region_curves |>
    filter(rl == "r", est_name != "Deconv-Tuned"),
  aes(time_value, est, color = est_name, alpha = est_name, linetype = est_name)
) +
  geom_line(linewidth = 0.7) +
  scale_alpha_manual(
    breaks = names(colours),
    values = c(1, rep(0.8, 10)),
    guide = "none"
  ) +
  scale_color_manual(
    breaks = c("Deconv-0", "Deconv-1", "Deconv-2"),
    values = colours
  ) +
  scale_linetype_manual(
    breaks = c("Deconv-0", "Deconv-1", "Deconv-2"),
    values = linetypes
  ) +
  theme_bw() +
  ylab("HFR") +
  xlab("") +
  scale_y_continuous(limits = ylim) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.margin = margin(t = -10)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

plot_grid(
  p1,
  p2,
  nrow = 1,
  rel_widths = c(0.5, 0.5),
  align = "h",
  labels = "AUTO"
)
ggsave(
  here::here("Figures", "BW", "Real_Time", "single_curve.pdf"),
  width = 8,
  height = 4
)

# tf_dates <- region_curves |> filter(est_name=="Deconv-0") |> pull(time_value)
gt_dates <- region_curves |>
  filter(est_name == "Ground Truth") |>
  pull(time_value) |>
  unique()
cr_dates <- region_curves |>
  filter(est_name == "Conv Ratio") |>
  pull(time_value)
cr_gt_dates <- intersect_dates(cr_dates, gt_dates)
region_gt <- region_curves |>
  filter(est_name == "Ground Truth", time_value %in% cr_gt_dates, rl == "l") |>
  pull(est)
region_tf0 <- region_curves |>
  filter(est_name == "Deconv-0", time_value %in% cr_gt_dates, rl == "r") |>
  pull(est)
region_tf1 <- region_curves |>
  filter(est_name == "Deconv-1", time_value %in% cr_gt_dates, rl == "r") |>
  pull(est)
region_tf2 <- region_curves |>
  filter(est_name == "Deconv-2", time_value %in% cr_gt_dates, rl == "r") |>
  pull(est)
region_tft <- region_curves |>
  filter(est_name == "Deconv-Tuned", time_value %in% cr_gt_dates) |>
  pull(est)
region_cr <- region_curves |>
  filter(est_name == "Conv Ratio", time_value %in% cr_gt_dates) |>
  pull(est)
mae(region_tf0, region_gt)
mae(region_tf1, region_gt)
mae(region_tf2, region_gt)
mae(region_tft, region_gt)
mae(region_cr, region_gt)

# ======================================================
# Figures for misspecified delay distributions
# ======================================================

# Retro misspecification
offset_map <- c(
  `1` = 0,
  `4` = -3,
  `5` = -2,
  `6` = -1,
  `7` = 1,
  `8` = 2,
  `9` = 3
)
misp_maes <- read_rds(here::here(
  "Data",
  "batch-results",
  "retro-maes-mis.rds"
)) |>
  mutate(
    est_name = case_when(
      is.na(order) ~ estimator,
      order == 0 ~ "Deconv-0",
      order == 1 ~ "Deconv-1",
      order == 2 ~ "Deconv-2"
    ),
    est_name = case_when(
      estimator == "lagged" ~ "Lagged Ratio",
      estimator == "conv" ~ "Conv Ratio",
      TRUE ~ est_name
    )
  )
tuned_misp_tf <- misp_maes |>
  # filter(algorithm == "trendfilter") |>
  filter(estimator == "trendfilter") |>
  summarize(
    mae = mae[cverr == min(cverr, na.rm = TRUE)],
    est_name = "Deconv-Tuned",
    .by = c(geo_value, mod_idx)
  )

conv_misp_mae <- misp_maes |>
  filter(estimator == "conv") |>
  select(conv_mae = mae, geo_value, mod_idx)

right_df <- bind_rows(
  misp_maes |> select(mae, geo_value, mod_idx, est_name),
  tuned_misp_tf
) |>
  left_join(conv_misp_mae, by = join_by(geo_value, mod_idx)) |>
  filter(mod_idx %in% as.integer(names(offset_map))) |>
  mutate(offset = offset_map[as.character(mod_idx)]) |>
  filter(est_name != "Deconv-Tuned")
left_df <- right_df |>
  mutate(value = (1 - mae / conv_mae) * 100) |>
  filter(!(est_name %in% c("Lagged Ratio", "Conv Ratio"))) |>
  select(value, est_name, offset)
right_df <- right_df |>
  rename(value = mae) |>
  select(value, est_name, offset)

retro_left <- left_df %>% mutate(source = "Retrospective")
retro_right <- right_df %>% mutate(source = "Retrospective")

# Real-time misspecification
rt_mis_maes <- read_rds(here::here("Data", "batch-results", "rt-maes-mis.rds"))

offset_map <- c(0, 0, 0, -5, -3, -1, 1, 3, 5)
names(offset_map) <- 1:9
offset_map <- offset_map[c(1, 4:9)] # ignore SD misspecifications

rt_mis_maes <- rt_mis_maes |>
  filter(mod_idx %in% as.numeric(names(offset_map))) |>
  mutate(
    offset = offset_map[as.character(mod_idx)],
    est_name = case_when(
      is.na(order) ~ estimator,
      order == 0 ~ "Deconv-0",
      order == 1 ~ "Deconv-1",
      order == 2 ~ "Deconv-2"
    ),
    est_name = case_when(
      est_name == "lagged" ~ "Lagged Ratio",
      est_name == "conv" ~ "Conv Ratio",
      est_name == "tf_tuned_order" ~ "Deconv-Tuned",
      TRUE ~ est_name
    )
  ) |>
  filter(est_name != "Deconv-Tuned")

conv_rt_mis_maes <- rt_mis_maes |>
  filter(est_name == "Conv Ratio") |>
  select(conv_mae = mae, geo_value, offset)

right_df <- rt_mis_maes |>
  select(mae, geo_value, offset, est_name) |>
  left_join(conv_rt_mis_maes, by = join_by(geo_value, offset))
left_df <- right_df |>
  mutate(value = (1 - mae / conv_mae) * 100) |>
  filter(!(est_name %in% c("Lagged Ratio", "Conv Ratio"))) |>
  select(value, est_name, offset)
right_df <- right_df |>
  rename(value = mae) |>
  select(value, est_name, offset)

rt_left <- left_df %>% mutate(source = "Real-Time")
rt_right <- right_df %>% mutate(source = "Real-Time")

# Combine retro + real-time
all_left <- bind_rows(retro_left, rt_left)
all_right <- bind_rows(retro_right, rt_right)

# Figure: MAE
all_right %>%
  mutate(source = factor(source, levels = c("Retrospective", "Real-Time"))) %>%
  group_by(source, offset, est_name) %>%
  summarise(
    m = mean(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE) / sqrt(51),
    .groups = "drop"
  ) %>%
  ggplot(aes(offset, m, color = est_name, linetype = est_name, shape = est_name)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.2) +
  facet_wrap(~source, nrow = 1, scales = "free_x") + # retro left, RT right
  scale_x_continuous(breaks = seq(-5, 5, 1), minor_breaks = seq(-5, 5, 1)) +
  ylab("MAE") +
  xlab("Mean Offset") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  ) +
  scale_color_manual(
    values = colours,
    guide = guide_legend(nrow = 1),
    breaks = methods
  ) +
  scale_linetype_manual(
    values = linetypes,
    guide = guide_legend(nrow = 1),
    breaks = methods
  ) +
  scale_shape_manual(
    values = shapes,
    guide = guide_legend(nrow = 1),
    breaks = methods
  )

ggsave(here::here("Figures", "BW", "misspecified_maes.pdf"), width = 8, height = 4)

# Figure: Percent improvements
all_left %>%
  mutate(source = factor(source, levels = c("Retrospective", "Real-Time"))) %>%
  group_by(source, offset, est_name) %>%
  summarise(
    m = mean(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE) / sqrt(51),
    .groups = "drop"
  ) %>%
  ggplot(aes(offset, m, color = est_name, linetype = est_name, shape = est_name)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.2) +
  facet_wrap(~source, nrow = 1, scales = "free_x") +
  scale_x_continuous(breaks = seq(-5, 5, 1), minor_breaks = seq(-5, 5, 1)) +
  ylab("Improv to Conv Ratio (%)") +
  xlab("Mean Offset") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  ) +
  scale_color_manual(values = colours, guide = guide_legend(nrow = 1)) +
  scale_linetype_manual(values = linetypes, guide = guide_legend(nrow = 1)) +
  scale_shape_manual(values = shapes, guide = guide_legend(nrow = 1))

ggsave(
  here::here("Figures", "BW", "misspecified_improvements.pdf"),
  width = 8,
  height = 4
)

# Figure - Real Data ----------------------------------------------------


source(here::here("Code", "funs.R"))
pa_curves <- read_rds(here::here("Data", "batch-results", "real-data-curves.rds"))

real_curves_retro <- pa_curves$retro_1se
# real_curves_rt <- pa_curves$rt_min_1se
real_curves_rt <- pa_curves$rt_1se_1se
library(forcats)
library(stringr)

plot_dat <- bind_rows(
  Retrospective = real_curves_retro,
  `Real-Time`   = real_curves_rt,
  .id = "type"
) |>
  transmute(
    type,
    time_value,
    est_name  = recode(
      estimator,
      conv   = "Conv Ratio",
      lagged = "Lagged Ratio",
      tf0    = "Deconv-0",
      tf1    = "Deconv-1",
      tf2    = "Deconv-2",
      .default = estimator
    ),
    est
  ) |>
  mutate(
    est_name = fct_relevel(
      est_name,
      "Conv Ratio", "Lagged Ratio", "Deconv-0", "Deconv-1", "Deconv-2"
    ),
    type = fct_relevel(type, "Retrospective")
  )

plot_dat |> group_by(type, est_name) |>
  summarize(min(time_value), max(time_value))

ylim <- range(plot_dat$est)

p1 <- ggplot(
  plot_dat |> drop_na() |>
    filter(type == "Retrospective",
           est_name %in% c("Conv Ratio", "Lagged Ratio", "Deconv-2")),
  aes(time_value, est, color = est_name, alpha = est_name, linetype = est_name)
) +
  geom_line(linewidth = 0.7) +
  scale_alpha_manual(
    breaks = names(colours),
    values = c(1, rep(0.8, 10)),
    guide = "none"
  ) +
  scale_color_manual(
    breaks = c("Conv Ratio", "Lagged Ratio", "Deconv-2"),
    values = colours
  ) +
  scale_linetype_manual(
    breaks = c("Conv Ratio", "Lagged Ratio", "Deconv-2"),
    values = linetypes
  ) +
  theme_bw() +
  ggtitle("Retrospective") +
  ylab("Estimated HFR") +
  xlab("") +
  scale_y_continuous(limits = ylim) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.margin = margin(t = -10)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

p2 <- ggplot(
  plot_dat |> drop_na() |>
    filter(type == "Real-Time",
           est_name %in% c("Conv Ratio", "Lagged Ratio", "Deconv-0")),
  aes(time_value, est, color = est_name, alpha = est_name, linetype = est_name)
) +
  geom_line(linewidth = 0.7) +
  scale_alpha_manual(
    breaks = names(colours),
    values = c(1, rep(0.8, 10)),
    guide = "none"
  ) +
  scale_color_manual(
    breaks = c("Conv Ratio", "Lagged Ratio", "Deconv-0"),
    values = colours
  ) +
  scale_linetype_manual(
    breaks = c("Conv Ratio", "Lagged Ratio", "Deconv-0"),
    values = linetypes
  ) +
  theme_bw() +
  ggtitle("Real-Time") +
  ylab("Estimated HFR") +
  xlab("") +
  scale_y_continuous(limits = ylim) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.margin = margin(t = -10)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

plot_grid(
  p1,
  p2,
  nrow = 1,
  rel_widths = c(0.5, 0.5),
  align = "h",
  labels = "AUTO"
)

ggsave(
  here::here("Figures", "BW", "real-data-pa.pdf"),
  width = 8,
  height = 4
)

# Appendix ----------------------------------------------------------------

# Figure, Misspecified delays -----------------------------------------------------

delay_means <- bind_rows(
  Retrospective = meta_data |>
    transmute(
      geo_value = tolower(geo_value),
      lag = retro_lag_argmax
    ),
  `Real-Time` = meta_data |>
    transmute(
      geo_value = tolower(geo_value),
      lag = rt_lag_argmax
    ),
  .id = "type"
) |>
  mutate(type = fct_relevel(factor(type), "Retrospective"))

# delay_means <- bind_rows(
#   Retrospective = c(
#     "US" = list(readRDS(here::here(
#       "Data",
#       "Simulated_Data",
#       'pb_sim_us.RData'
#     ))),
#     readRDS(here::here("Data", "Simulated_Data", 'pb_sim_state.RData'))
#   ) |>
#     map_dbl("Lag") |>
#     enframe(name = "geo_value", value = "lag") |>
#     mutate(geo_value = tolower(geo_value)),
#   `Real-Time` = read_rds(here::here(
#     "Data",
#     "National_Data",
#     "archive_summary_stats.rds"
#   )) |>
#     select(geo_value, lag),
#   .id = "type"
# ) |>
#   mutate(type = fct_relevel(factor(type), "Retrospective"))

locs <- c("ca", "la", "wy")
delay_means <- filter(delay_means, geo_value %in% locs)

true_delays <- delay_means |>
  rowwise() |>
  mutate(
    zz = list(tibble(
      x = 0:60,
      distn = make_delay_distr(lag, lag * 0.9, 60)
    ))
  ) |>
  ungroup() |>
  select(type, geo_value, zz) |>
  unnest(zz)

misp_delays <- delay_means |>
  mutate(
    # mean_modifier = list(c(0, 0, -5, -3, -1, 1, 3, 5)),
    # sd_mult = list(c(1, .8, rep(.9, 6)))
    mean_modifier = list(c(-5, -3, -1, 1, 3, 5)),
    sd_mult = list(rep(.9, 6))
  ) |>
  unnest_longer(c(mean_modifier, sd_mult)) |>
  rowwise() |>
  mutate(
    zz = list(tibble(
      x = 0:60,
      distn = make_delay_distr(
        pmax(lag + mean_modifier, 1),
        pmax(lag + mean_modifier, 1) * sd_mult,
        60L
      )
    ))
  ) |>
  select(type, geo_value, zz) |>
  ungroup() |>
  mutate(id = row_number(), .by = c(geo_value, type)) |>
  unnest(zz)

ggplot(misp_delays, aes(x, distn)) +
  geom_line(aes(color = factor(id), group = factor(id))) +
  geom_line(data = true_delays, color = "black", linewidth = 1.2) +
  facet_grid(
    type ~ geo_value,
    scales = "free_y",
    labeller = labeller(geo_value = toupper),
    switch = "y"
  ) +
  theme_bw() +
  scale_color_grey(start = 0.75, end = 0.1) +
  scale_y_continuous(expand = expansion(c(0, 0.05)), name = "") +
  scale_x_continuous(expand = expansion(0), name = "", breaks = c(0, 20, 40)) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(hjust = 0)
  )

ggsave(
  here::here("Figures", "BW", "misspecified_delay_distrs.pdf"),
  width = 8,
  height = 4
)

# Figure - Retrospective HFR curves ------------------------------------

# Constant-order trend filtering
all_retro_curves <- retro_curves_og |>
  filter(repl==1) |>
  mutate(
    est_name = case_when(
      estimator == "conv" ~ "Conv Ratio",
      estimator == "lagged" ~ "Lagged Ratio",
      order == 0 ~ "Deconv-0",
      order == 1 ~ "Deconv-1",
      order == 2 ~ "Deconv-2"
    )
  ) |>
  rename(geo_value = region) |>
  unnest(result) |>
  select(time_value, est, est_name, geo_value) |>
  bind_rows(
    gt_hfrs |>
      mutate(est_name = "Ground Truth") |>
      rename(est = hfr)
  ) |>
  mutate(
    est_name = fct_relevel(
      as.factor(est_name),
      "Ground Truth",
      "Conv Ratio",
      "Lagged Ratio",
      "Deconv-0",
      "Deconv-1",
      "Deconv-2"
    ),
    ww = ifelse(est_name == "Ground Truth", "a", "b")
  )

# Only plot during estimated period (excluding burn-in for TF)
conv_dates <- all_retro_curves |>
  filter(est_name == "Conv Ratio") |>
  group_by(geo_value) |>
  summarise(
    first_date = min(time_value),
    last_date = max(time_value),
    .groups = "drop"
  )

all_retro_curves <- all_retro_curves |>
  left_join(conv_dates, by = "geo_value") |>
  filter(time_value >= first_date, time_value <= last_date) |>
  select(-c(first_date, last_date))

all_retro_curves <- all_retro_curves |>
  mutate(est_name = str_replace(est_name, "^TF", "Deconv"))

ggplot(
  all_retro_curves |> 
    filter(geo_value != "us") |>
    filter(est_name %in% c(
      "Conv Ratio",
      "Lagged Ratio",
      #"Deconv-0",
      #"Deconv-1",
      "Deconv-2"
      )) |>
    drop_na(est) |>
    group_by(geo_value, est_name) |>
    slice(seq(1, n(), by = 7)) |>
    ungroup(),
  aes(time_value, est, color = est_name, alpha = est_name, linetype = est_name)
) +
  geom_line() +
  facet_wrap(~geo_value, ncol = 5, labeller = labeller(geo_value = toupper)) +
  coord_cartesian(ylim = c(.05, .45)) +
  scale_alpha_manual(
    breaks = names(colours),
    values = c(1, rep(0.8, 10)),
    guide = "none"
  ) +
  scale_color_manual(
    breaks = c(
      "Conv Ratio",
      "Lagged Ratio",
      #"Deconv-0",
      #"Deconv-1",
      "Deconv-2"
    ),
    values = colours
  ) +
  scale_linetype_manual(
    breaks = c("Conv Ratio", "Lagged Ratio", "Deconv-2"),
    values = linetypes
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +
  ylab("HFR") +
  xlab("") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.margin = margin(t = -10)
  )

ggsave(
  here::here("Figures", "BW", "Retrospective", "all_curves.pdf"),
  width = 11,
  height = 15,
  dpi = 300
)

# Figure - Real Time HFR curves ------------------------------------
# 0-order trend filtering
all_rt_curves <- rt_curves |>
  rename(geo_value = region) |>
  filter(estimator != "tf_tuned_order") |>
  mutate(
    est_name = case_when(
      estimator == "conv" ~ "Conv Ratio",
      estimator == "lagged" ~ "Lagged Ratio",
      estimator == "trendfilter" & order == 0 ~ "Deconv-0",
      estimator == "trendfilter" & order == 1 ~ "Deconv-1",
      estimator == "trendfilter" & order == 2 ~ "Deconv-2",
    )
  ) |>
  bind_rows(
    gt_hfrs |>
      mutate(est_name = "Ground Truth") |>
      rename(est = hfr)
  ) |>
  mutate(
    est_name = fct_relevel(
      as.factor(est_name),
      "Ground Truth",
      "Conv Ratio",
      "Lagged Ratio",
      "Deconv-0",
      "Deconv-1",
      "Deconv-2"
    ),
    ww = ifelse(est_name == "Ground Truth", "a", "b")
  )
# 
# # Subset to proper region
# fourth_conv_dates <- all_rt_curves |>
#   filter(est_name == "Conv Ratio") |>
#   group_by(geo_value) |>
#   summarise(first_date = min(time_value))[4], .groups = "drop")
# 
# last_gt_dates <- all_rt_curves |>
#   filter(est_name == "Ground Truth") |>
#   group_by(geo_value) |>
#   summarise(last_date = max(time_value), .groups = "drop")
# 
# all_rt_curves <- all_rt_curves |>
#   left_join(fourth_conv_dates, by = "geo_value") |>
#   filter(time_value >= fourth_date) |>
#   select(-fourth_date) |>
#   left_join(last_gt_dates, by = "geo_value") |>
#   filter(time_value <= last_date) |>
#   select(-last_date)

# Only plot during estimated period (excluding burn-in for TF)
conv_dates <- all_rt_curves |>
  filter(est_name == "Conv Ratio") |>
  group_by(geo_value) |>
  summarise(
    first_date = min(time_value),
    last_date = max(time_value),
    .groups = "drop"
  )

all_rt_curves <- all_rt_curves |>
  left_join(conv_dates, by = "geo_value") |>
  filter(time_value >= first_date, time_value <= last_date) |>
  select(-c(first_date, last_date))

all_rt_curves <- all_rt_curves |>
  mutate(est_name = str_replace(est_name, "^TF", "Deconv"))

ggplot(
  all_retro_curves |> 
    filter(geo_value != "us") |>
    filter(est_name %in% c(
      "Conv Ratio",
      "Lagged Ratio",
      "Deconv-0"
      #"Deconv-1",
      #"Deconv-2")
      )) |>
    drop_na(est),
  aes(time_value, est, color = est_name, alpha = est_name, linetype = est_name)
) +
  geom_line() +
  facet_wrap(~geo_value, ncol = 5, labeller = labeller(geo_value = toupper)) +
  coord_cartesian(ylim = c(.05, .45)) +
  scale_alpha_manual(
    breaks = names(colours),
    values = c(1, rep(0.8, 10)),
    guide = "none"
  ) +
  scale_color_manual(
    breaks = c(
      "Conv Ratio",
      "Lagged Ratio",
      "Deconv-0"
      #"Deconv-1",
      #"Deconv-2"
    ),
    values = colours
  ) +
  scale_linetype_manual(
    breaks = c("Conv Ratio", "Lagged Ratio", "Deconv-0"),
    values = linetypes
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +
  ylab("HFR") +
  xlab("") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.margin = margin(t = -10)
  )

ggsave(
  here::here("Figures", "BW", "Real_Time", "all_curves.pdf"),
  width = 11,
  height = 15,
  dpi=300
)



# Figure - Tail regularization comparison ---------------------------------

tail_reg_curves <- read_rds(here::here("Data", "batch-results", "tail-reg-curves.rds"))
est_date <- sort(unique(
  read_rds(here::here("Data", "raw-data", "full-archive.rds"))$version
))[90L]

tail_reg_curves %>%
  filter(time_value >= min(tail_reg_curves$time_value) + 4 * 60) %>%
  filter(order != "Order 0 (Constant)") %>%
  ggplot(aes(time_value, value, color = gamma, group = gamma)) +
    geom_line() +
    # facet_wrap(~order, scales = "free_y") +
    facet_wrap(~order) +
  labs(
    title = expression("Deconvolution HFR estimates by "* gamma),
    subtitle = paste0("Estimating at ", est_date, ". Fixing lambda."),
    x = "Date", y = "HFR"
  ) +
  theme_bw() +
  scale_color_gradient(
    trans = "log10",
    low = "grey85",
    high = "grey10",
    name = expression(gamma)
  ) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    plot.title = element_text(size = 16),#, face = "bold"
    plot.subtitle = element_text(size = 13),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11),
    strip.text = element_text(size = 12)
  )
ggsave(here::here("Figures", "BW", "tail-regularization-comparison.pdf"),
       width = 10, height = 5)



# Figure - Runtime comparison ---------------------------------

timing_df <- read_rds(here::here("Data", "benchmark-timings.rds"))

timing_df %>%
  mutate(expression = recode(expression,
                             "tf0"    = "Deconv-0",
                             "tf1"    = "Deconv-1",
                             "tf2"    = "Deconv-2",
                             "conv"   = "Conv Ratio",
                             "lagged" = "Lagged Ratio"
  )) %>%
  mutate(expression = factor(expression, levels = c("Conv Ratio", "Lagged Ratio", "Deconv-0", "Deconv-1", "Deconv-2"))) %>%
  group_by(expression, n) %>%
  summarise(
    mean_time = mean(as.numeric(time)),
    se_time   = sd(as.numeric(time)) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(
    x = n, y = mean_time,
    color = expression, group = expression,
    linetype = expression, shape = expression
  )) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_time - se_time, ymax = mean_time + se_time), width = 0.05) +
  scale_x_log10(breaks = unique(timing_df$n)) +
  scale_color_manual(values = colours) +
  scale_linetype_manual(values = linetypes) +
  scale_shape_manual(values = shapes) +
  labs(
    title = "Runtime comparison for retrospective estimation",
    # x     = "Number of hospitalization dates",
    x     = "Number of hospitalization timesteps",
    y     = "Time (seconds)",
    color = "Method",
    linetype = "Method",
    shape = "Method"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  
  theme(
    # legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    plot.title = element_text(size = 16),#, face = "bold"
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11),
    legend.text = element_text(size = 11),
    strip.text = element_text(size = 12)
  )

ggsave(here::here("Figures", "BW", "runtime-comparison.pdf"),
       width = 8, height = 4)


# Time-varying delay experiment (Appendix) --------------------------------
# Black-and-white versions. These read committed results (MAE summary produced
# by Code/Analysis/maes-time-var.R); each block runs in local() so its
# temporaries do not leak into the global environment.
tvar_out <- here::here("Figures", "BW", "time-varying")

# Figure - Time-varying MAE summary (greyscale fills)
local({
  summary_tab <- read_rds(here::here(
    "Data", "batch-results", "time-var-mae-summary.rds"
  ))
  method_levels <- c("tf-o0", "tf-o1", "tf-o2", "conv", "lagged")
  method_labels <- c("Deconv-0", "Deconv-1", "Deconv-2",
                     "Conv. ratio", "Lagged ratio")
  state_levels <- c("ca", "la", "pa", "ut")
  state_labels <- c("CA", "LA", "PA", "UT")
  plot_df <- summary_tab |>
    mutate(
      method = factor(method, levels = method_levels, labels = method_labels),
      geo_value = factor(geo_value, levels = state_levels, labels = state_labels),
      problem = factor(problem, levels = c("retro", "rt"),
                       labels = c("Retrospective", "Real-time")),
      spec = factor(spec, levels = c("ws", "mis"),
                    labels = c("Well-specified delay", "Misspecified delay"))
    )
  tvar_df <- plot_df |>
    filter(experiment == "tvar") |>
    mutate(se = sd_mae / sqrt(n), lo = mean_mae - se, hi = mean_mae + se)
  p <- ggplot(tvar_df, aes(x = geo_value, y = mean_mae, fill = method)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.78,
             color = "black", linewidth = 0.2) +
    geom_errorbar(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.85),
                  width = 0.35, linewidth = 0.3, color = "grey25") +
    facet_grid(problem ~ spec, scales = "free_y") +
    scale_fill_manual(
      name = NULL,
      values = c(
        "Deconv-0" = "grey0", "Deconv-1" = "grey30", "Deconv-2" = "grey50",
        "Conv. ratio" = "grey70", "Lagged ratio" = "grey88"
      )
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Mean absolute error",
         caption = "Mean MAE +/- 1 SE across 10 simulation replicates.") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          plot.caption = element_text(hjust = 0, size = 8),
          panel.grid.major.x = element_blank())
  ggsave(file.path(tvar_out, "mae-summary.pdf"), p, width = 8, height = 4.5)
})

# Figure - Time-varying curves (retro, well-specified; greyscale + linetype)
local({
  res <- read_rds(here::here("Data", "batch-results", "time-var-delay.rds"))
  gt <- read_rds(here::here("Data", "raw-data", "daily-hfrs.rds")) |>
    filter(geo_value %in% c("ca", "la", "pa", "ut")) |>
    rename(truth = hfr)
  res <- res |>
    mutate(method = case_when(
      algorithm %in% c("tf-rt", "tf-retro") ~ paste0("tf-o", order),
      estimator == "conv" ~ "conv",
      estimator == "lagged" ~ "lagged",
      TRUE ~ paste(algorithm, estimator, sep = "-")
    ))
  med_curves <- res |>
    select(problem, region, method, spec, repl, res) |>
    unnest(res) |>
    rename(geo_value = region) |>
    summarise(est = median(est, na.rm = TRUE),
              .by = c(problem, geo_value, method, spec, time_value))
  method_levels <- c("tf-o0", "tf-o1", "tf-o2", "conv")
  method_labels <- c("Deconv-0", "Deconv-1", "Deconv-2", "Conv. ratio")
  state_levels <- c("ca", "la", "pa", "ut")
  state_labels <- c("CA", "LA", "PA", "UT")
  method_greys <- c("Deconv-0" = "grey0", "Deconv-1" = "grey35",
                    "Deconv-2" = "grey55", "Conv. ratio" = "grey45")
  method_ltys <- c("Deconv-0" = "solid", "Deconv-1" = "dashed",
                   "Deconv-2" = "dotdash", "Conv. ratio" = "dotted")
  date_start <- as.Date("2021-01-01")
  date_end <- as.Date("2022-04-01")
  d <- med_curves |>
    filter(problem == "retro", spec == "ws",
           time_value >= date_start, time_value <= date_end,
           method %in% method_levels) |>
    mutate(method = factor(method, levels = method_levels, labels = method_labels),
           geo_value = factor(geo_value, levels = state_levels, labels = state_labels))
  gt_d <- gt |>
    filter(time_value >= date_start, time_value <= date_end) |>
    mutate(geo_value = factor(geo_value, levels = state_levels, labels = state_labels))
  p <- ggplot(d, aes(time_value, est, color = method, linetype = method)) +
    geom_line(data = gt_d, aes(time_value, truth),
              color = "black", linewidth = 0.7, linetype = "solid",
              inherit.aes = FALSE) +
    geom_line(linewidth = 0.55) +
    facet_wrap(~ geo_value, ncol = 2, scales = "free_y") +
    scale_color_manual(values = method_greys, name = NULL) +
    scale_linetype_manual(values = method_ltys, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    guides(color = guide_legend(override.aes = list(linewidth = 1.0))) +
    labs(x = NULL, y = "HFR",
         title = "Retrospective estimates under well-specified time-varying delay",
         subtitle = "Median across 10 simulation replicates. Black line: ground truth.") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 14),
          legend.key.width = unit(0.9, "cm"),
          legend.key.height = unit(0.55, "cm"),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 10),
          panel.grid.minor = element_blank())
  ggsave(file.path(tvar_out, "curves-retro-ws.pdf"), p, width = 9, height = 4.5)
})

# Figure - Time-varying delay illustration (greyscale + linetype)
local({
  delay_means <- tribble(
    ~variant,   ~retro, ~rt,
    "Original",  8,     18,
    "Alpha",     5,      6,
    "Delta",    13,     22,
    "Omicron",   9,     21
  )
  state_levels <- c("ca", "la", "pa", "ut")
  state_labels <- c("CA", "LA", "PA", "UT")
  variant_levels <- c("Original", "Alpha", "Delta", "Omicron")
  variant_greys <- c("Original" = "grey80", "Alpha" = "grey55",
                     "Delta" = "grey30", "Omicron" = "grey0")
  variant_ltys <- c("Original" = "solid", "Alpha" = "dashed",
                    "Delta" = "dotdash", "Omicron" = "dotted")
  date_start <- as.Date("2021-01-01")
  date_end   <- as.Date("2022-04-01")
  dd <- 60L
  pmfs <- delay_means |>
    rowwise() |>
    mutate(retro_pmf = list(make_delay_distr(retro, 0.9 * retro, dd)),
           rt_pmf    = list(make_delay_distr(rt,    0.9 * rt,    dd))) |>
    ungroup() |>
    mutate(day = list(0:dd)) |>
    select(variant, retro_pmf, rt_pmf, day) |>
    unnest(c(retro_pmf, rt_pmf, day)) |>
    pivot_longer(c(retro_pmf, rt_pmf), names_to = "problem", values_to = "prob") |>
    mutate(problem = recode(problem, retro_pmf = "Retrospective", rt_pmf = "Real-time"),
           problem = factor(problem, levels = c("Retrospective", "Real-time")),
           variant = factor(variant, levels = variant_levels))
  props <- read_rds(here::here("Data", "raw-data", "daily-variant-props.rds")) |>
    filter(geo_value %in% state_levels,
           time_value >= date_start, time_value <= date_end) |>
    pivot_longer(Original:Omicron, names_to = "variant", values_to = "prop") |>
    mutate(variant = factor(variant, levels = variant_levels),
           geo_value = factor(geo_value, levels = state_levels, labels = state_labels))
  plot_max_day <- 45
  p_top <- ggplot(pmfs |> filter(day <= plot_max_day),
                  aes(day, prob, color = variant, linetype = variant)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ problem, nrow = 1) +
    scale_color_manual(values = variant_greys, name = NULL) +
    scale_linetype_manual(values = variant_ltys, name = NULL) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Delay (days from hospitalization to death)", y = "Probability",
         title = "(a) Per-variant delay distributions (shared across the four states)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.text = element_text(size = 12),
          legend.key.width = unit(1.0, "cm"), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 12))
  p_bot <- ggplot(props, aes(time_value, prop, fill = variant)) +
    geom_area(position = "stack", color = "grey30", linewidth = 0.15) +
    facet_wrap(~ geo_value, nrow = 1) +
    scale_fill_manual(values = variant_greys, name = NULL, guide = "none") +
    scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(accuracy = 1)) +
    scale_x_date(expand = c(0, 0),
                 breaks = as.Date(c("2021-01-01", "2021-07-01", "2022-01-01")),
                 date_labels = "%Y-%m", limits = c(date_start, date_end)) +
    labs(x = NULL, y = "Variant share",
         title = "(b) Daily variant composition per state") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 12))
  p <- p_top / p_bot + plot_layout(heights = c(1.1, 1))
  ggsave(file.path(tvar_out, "delay-illustration.pdf"), p, width = 10, height = 5)
})
