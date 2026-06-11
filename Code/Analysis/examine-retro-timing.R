library(tidyverse)
library(data.table)
library(bench)

source(here::here("Code", "funs.R"))
source(here::here("Code", "helper_functions.R"))
source(here::here("Code", "estimate.R"))
source(here::here("Code", "batch-funs.R"))


# Global parameters -------------------------------------------------------

d <- 60L
ws <- c(1, 7, 15, 21, 29)
ord <- 0:2
regions <- c("pa", "ut", "ks", "us") # each has 964 dates
ns <- c(150, 200, 300, 500, 750, 964)

# Grab data ---------------------------------------------------------------
data_finalized <- read_rds(here::here("Data", "raw-data", "finalized.rds")) |>
  filter(geo_value %in% regions)
meta_data <- read_rds(here::here("Data", "raw-data", "all-meta-data.rds")) |>
  filter(geo_value %in% regions)

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

create_data_problem <- function(region, n) {
  local_data <- data_finalized |>
    filter(geo_value == region) |>
    arrange(time_value) |>
    slice_head(n = n)
  delay_mean <- meta_data |>
    filter(geo_value == region) |>
    pull(retro_lag_argmax)
  delay_sd <- delay_mean * 0.9
  delay_distn <- make_delay_distr(delay_mean, delay_sd, d)
  ll <- process_tibble_to_named_vctrs(local_data)
  ndeath <- length(ll$dd)
  nhosp <- length(ll$hh)

  Z <- make_Z_mat_simple(ndeath, ll$hh, delay_distn)
  Y <- ll$dd
  maxmax <- max(max(abs(Z)), max(abs(Y)))
  Z <- Z / maxmax
  Y <- Y / maxmax
  c(
    ll,
    Z = list(Z),
    Y = list(Y),
    delay_distn = list(delay_distn),
    delay_mean = list(delay_mean),
    lambda = 100,
    w = 14
  )
}

tf_single_lambda <- function(Z, Y, delay_distr, lambda, order = 2) {
  out <- tryCatch(
    {
      compute_tf_hfrs(
        Z,
        Y,
        lambda = lambda,
        order = order,
        real_time = FALSE,
        d = d,
        delay_distr_T = delay_distr
      )
    },
    error = function(cond) {
      return(rep(NA, ncol(Z)))
    }
  )
  if (is.null(out)) {
    out <- rep(NA, ncol(Z))
  }
  out
}


timings <- bench::press(
  region = regions,
  n = ns,
  {
    dat <- create_data_problem(region, n)
    bench::mark(
      check = FALSE,
      iterations = 20,
      lagged = compute_lagged_hfrs_fast(
        dat$hh,
        dat$dd,
        lag = dat$delay_mean,
        w = dat$w,
        real_time = FALSE
      ),
      conv = compute_conv_hfrs_retro(dat$hh, dat$dd, dat$delay_distn),
      tf0 = tf_single_lambda(
        dat$Z,
        dat$Y,
        dat$delay_distn,
        dat$lambda,
        order = 0
      ),
      tf1 = tf_single_lambda(
        dat$Z,
        dat$Y,
        dat$delay_distn,
        dat$lambda,
        order = 1
      ),
      tf2 = tf_single_lambda(
        dat$Z,
        dat$Y,
        dat$delay_distn,
        dat$lambda,
        order = 2
      )
    )
  }
)

timings <- timings |> select(expression, region, n, time)
write_rds(timings, here::here("Data", "benchmark-timings.rds"))

df <- read_rds(here::here("Data", "benchmark-timings.rds"))

df %>%
  mutate(
    expression = recode(
      expression,
      "tf0" = "Deconv-0",
      "tf1" = "Deconv-1",
      "tf2" = "Deconv-2",
      "conv" = "Conv Ratio",
      "lagged" = "Lagged Ratio"
    )
  ) %>%
  group_by(expression, n) %>%
  summarise(
    mean_time = mean(as.numeric(time)),
    se_time = sd(as.numeric(time)) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = n, y = mean_time, color = expression, group = expression)) +
  geom_line() +
  geom_point() +
  geom_errorbar(
    aes(ymin = mean_time - se_time, ymax = mean_time + se_time),
    width = 0.05
  ) +
  scale_x_log10(breaks = unique(df$n)) +
  # scale_color_viridis_d() +
  scale_color_manual(
    values = c(
      "Conv Ratio" = "#7B2D8B", # purple
      "Lagged Ratio" = "#2166AC", # blue
      "Deconv-0" = "#1B7837", # dark green
      "Deconv-1" = "#5AAE61", # medium green
      "Deconv-2" = "#F4A736" # warm amber (instead of yellow)
    )
  ) +
  labs(
    title = "Runtime comparison",
    x = "Number of hospitalization dates",
    y = "Time (seconds)",
    color = "Method"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
