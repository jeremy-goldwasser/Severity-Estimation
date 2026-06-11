library(tidyverse)
library(epiprocess)
library(data.table)

source(here::here("Code", "funs.R"))
source(here::here("Code", "helper_functions.R"))
source(here::here("Code", "estimate.R"))
source(here::here("Code", "batch-funs.R"))


# Global parameters -------------------------------------------------------
d <- 60L
ord <- 2
region <- "pa"
est_date_loc <- 90L # which of the ~120 RT dates to use

# Grab data ---------------------------------------------------------------
data_archive <- read_rds(here::here("Data", "raw-data", "full-archive.rds")) |>
  as.data.frame() |>
  filter(geo_value == region) |>
  as_epi_archive()
meta_data <- read_rds(here::here("Data", "raw-data", "all-meta-data.rds")) |>
  filter(geo_value == region)
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

local_data <- data_archive_fixed |> epix_as_of(versions[est_date_loc])
delay_mean <- meta_data |> pull(rt_lag_1se)
delay_sd <- delay_mean * 0.9
delay_distn <- make_delay_distr(delay_mean, delay_sd, d)

tf_gamma_curves <- function(
  hosps,
  deaths,
  delay_distr,
  lambda_idx = 10,
  n_gammas = 20,
  order = 2
) {
  gammas = geom_spacing(1e0, 1e5, n_gammas)
  ndeath <- length(deaths)
  nhosp <- length(hosps)
  d <- length(delay_distr) - 1
  stopifnot(ndeath + d == nhosp)

  Z <- make_Z_mat_simple(ndeath, hosps, delay_distr)
  Y <- deaths
  maxmax <- max(max(abs(Z)), max(abs(Y)))
  Z <- Z / maxmax
  Y <- Y / maxmax

  lambda_max <- compute_lambda_max(Z, Y, order, family = "poissonLI")
  lambdas <- geom_spacing(lambda_max * 1e-6, lambda_max, 30)
  lambda <- lambdas[lambda_idx]

  full_fit <- map(gammas, \(gamma) {
    out <- tryCatch(
      {
        compute_tf_hfrs(
          Z,
          Y,
          lambda = lambda,
          order = order,
          real_time = TRUE,
          gamma = gamma,
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
    return(out)
  })
  full_fit <- map2(full_fit, gammas, \(x, y) {
    as_tibble(x) |> mutate(gamma = y, time_value = as.Date(names(hosps)))
  })
  full_fit |> list_rbind()
}

# Illustrate results ------------------------------------------------------
cache_file <- here::here("Data", "batch-results", "tail-reg-curves.rds")
if (file.exists(cache_file)) {
  curves <- read_rds(cache_file)
} else {
  set.seed(12345)
  ll <- process_tibble_to_named_vctrs(local_data)
  # curves_ord0 <- with(ll, tf_gamma_curves(hh, dd, delay_distn, lambda_idx = 23, order = 0)) %>%
  #   mutate(order = "Order 0 (Constant)")
  curves_ord1 <- with(
    ll,
    tf_gamma_curves(hh, dd, delay_distn, lambda_idx = 20, order = 1)
  ) %>%
    mutate(order = "Order 1 (Linear)")
  curves_ord2 <- with(
    ll,
    tf_gamma_curves(hh, dd, delay_distn, lambda_idx = 18, order = 2)
  ) %>%
    mutate(order = "Order 2 (Quadratic)")
  # curves <- bind_rows(curves_ord0, curves_ord1, curves_ord2) %>%
  #   mutate(order = factor(order, levels = c("Order 0 (Constant)", "Order 1 (Linear)", "Order 2 (Quadratic)")))
  curves <- bind_rows(curves_ord1, curves_ord2) %>%
    mutate(
      order = factor(
        order,
        levels = c("Order 1 (Linear)", "Order 2 (Quadratic)")
      )
    )
  write_rds(curves, cache_file)
}

est_date <- versions[est_date_loc]
curves %>%
  filter(time_value >= min(curves$time_value) + 4 * d) %>%
  filter(order != "Order 0 (Constant)") %>%
  ggplot(aes(time_value, value, color = gamma, group = gamma)) +
  geom_line() +
  facet_wrap(~order) +
  labs(
    title = expression("Deconvolution HFR estimates by " * gamma),
    subtitle = paste0("Estimating at ", est_date, ". Fixing lambda."),
    x = "Date",
    y = "HFR"
  ) +
  theme_bw() +
  scale_color_viridis_c(trans = "log10", name = expression(gamma)) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    plot.title = element_text(size = 18), #, face = "bold"
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11),
    strip.text = element_text(size = 12)
  )
ggsave(
  here::here("Figures", "tail-regularization-comparison.pdf"),
  width = 10,
  height = 5
)
