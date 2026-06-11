library(epidatr)
library(epiprocess)
library(tidyverse)
library(RcppRoll)
library(zoo)
# library(cli)

source(here::here("Code", "funs.R"))
d <- 60L
excl_geos <- c("as", "dc", "gu", "mp", "pr", "vi")
`%nin%` <- function(x, table) !(x %in% table)

first_date <- ymd("2020-03-01")
first_hosp_date <- ymd("2020-07-14")
first_hfr_date <- ymd("2020-08-08")
last_jhu <- ymd("2023-03-09") # case/death reporting ends
first_as_of <- ymd("2020-11-22") # first hhs issue in the api
first_as_of <- max(first_as_of, first_hosp_date + d)
est_rt_dates <- seq(first_as_of, last_jhu, by = 7)


# Download Real Time Data -----------------------------------------------------

deaths_archive <- list()
hosp_archive <- list()

for (i in seq_along(est_rt_dates)) {
  as_of <- est_rt_dates[i]
  cli::cli_inform("Downloading {as_of}, {i} out of {length(est_rt_dates)}.")
  jhu_st <- pub_covidcast(
    source = "jhu-csse",
    signals = "deaths_incidence_num",
    geo_type = "state",
    time_type = "day",
    as_of = as_of
  ) |>
    select(time_value, geo_value, deaths = value)
  jhu_nat <- pub_covidcast(
    source = "jhu-csse",
    signals = "deaths_incidence_num",
    geo_type = "nation",
    time_type = "day",
    as_of = as_of
  ) |>
    select(time_value, geo_value, deaths = value)
  jhu <- bind_rows(jhu_st, jhu_nat) |>
    filter(geo_value %nin% excl_geos) |>
    arrange(geo_value, time_value) |>
    complete(
      time_value = full_seq(c(time_value, as_of), period = 1),
      geo_value
    ) |>
    group_by(geo_value) |>
    filter(time_value >= ymd(first_hosp_date)) |>
    mutate(deaths = fix_right_tail(deaths)) |>
    ungroup()
  deaths_archive[[i]] <- jhu

  hhs_st <- pub_covidcast(
    source = "hhs",
    signals = "confirmed_admissions_covid_1d",
    geo_type = "state",
    time_type = "day",
    as_of = as_of
  ) |>
    select(time_value, geo_value, hosps = value)
  hhs_nat <- pub_covidcast(
    source = "hhs",
    signals = "confirmed_admissions_covid_1d",
    geo_type = "nation",
    time_type = "day",
    as_of = as_of
  ) |>
    select(time_value, geo_value, hosps = value)
  hhs <- bind_rows(hhs_st, hhs_nat) |>
    filter(geo_value %nin% excl_geos) |>
    arrange(geo_value, time_value) |>
    complete(
      time_value = full_seq(c(time_value, as_of), period = 1),
      geo_value
    ) |>
    group_by(geo_value) |>
    mutate(hosps = replace_na(hosps, 0)) |>
    filter(time_value >= ymd(first_hosp_date)) |>
    mutate(hosps = fix_right_tail(hosps)) |>
    ungroup()
  hosp_archive[[i]] <- hhs
}

deaths_archive <- deaths_archive |>
  set_names(est_rt_dates) |>
  list_rbind(names_to = "version")
hosp_archive <- hosp_archive |>
  set_names(est_rt_dates) |>
  list_rbind(names_to = "version")
full_archive <- full_join(
  deaths_archive,
  hosp_archive,
  by = join_by(time_value, version, geo_value)
) |>
  mutate(version = ymd(version)) |>
  arrange(version, geo_value, time_value) |>
  as_epi_archive()

write_rds(
  full_archive$DT,
  here::here("Data", "raw-data", "full-archive.rds")
)

# Download retrospective data -------------------------------------------------

nchs_st <- pub_covidcast(
  source = "nchs-mortality",
  signals = "deaths_covid_incidence_num",
  geo_type = "state",
  time_type = "week"
) |>
  select(time_value, geo_value, deaths = value)
nchs_nat <- pub_covidcast(
  source = "nchs-mortality",
  signals = "deaths_covid_incidence_num",
  geo_type = "nation",
  time_type = "week"
) |>
  select(time_value, geo_value, deaths = value)
nchs <- bind_rows(nchs_st, nchs_nat) |>
  filter(geo_value %nin% excl_geos) |>
  arrange(geo_value, time_value) |>
  # Deaths counted from Sunday to the following Saturday
  # Shift to use end-of-week: Saturday total of past week
  mutate(time_value = time_value + 6) |>
  filter(time_value >= ymd(first_hosp_date), time_value <= last_jhu) |>
  group_by(geo_value) |>
  # Interpolate internal missing values linearly.
  # Missingness is due to data suppression of low counts (1-9)
  mutate(deaths = na.approx(deaths, na.rm = FALSE)) |>
  # remove leading/trailing
  filter(!is.na(deaths)) |>
  # Now create daily: divide weekly totals by 7 and fill over the week
  complete(
    time_value = full_seq(c(time_value, min(time_value) - 6), period = 1)
  ) |>
  mutate(deaths = deaths / 7) |>
  fill(deaths, .direction = "up") |>
  ungroup()

hhs_st <- pub_covidcast(
  source = "hhs",
  signals = "confirmed_admissions_covid_1d",
  geo_type = "state",
  time_type = "day"
) |>
  select(time_value, geo_value, hosps = value)
hhs_nat <- pub_covidcast(
  source = "hhs",
  signals = "confirmed_admissions_covid_1d",
  geo_type = "nation",
  time_type = "day"
) |>
  select(time_value, geo_value, hosps = value)
hhs <- bind_rows(hhs_st, hhs_nat) |>
  filter(geo_value %nin% excl_geos) |>
  filter(time_value >= ymd(first_hosp_date), time_value <= last_jhu) |>
  arrange(geo_value, time_value) |>
  group_by(geo_value) |>
  mutate(hosps = replace_na(hosps, 0)) |>
  # only in beginning of ND, two dates in NE
  ungroup()

full_retro <- inner_join(
  nchs,
  hhs,
  by = join_by(time_value, geo_value)
) |>
  arrange(geo_value, time_value) |>
  relocate(geo_value, time_value)

write_rds(
  full_retro,
  here::here("Data", "raw-data", "finalized.rds")
)


# Create Naive HFRs -----------------------------------------------------------
vocs <- c("Original", "Alpha", "Delta", "Omicron")

us_variants <- read_csv(here::here(
  "Data",
  "Variants",
  "seq_df_us_biweekly.csv"
)) |>
  mutate(Other = case_when(Date >= ymd("2022-01-01") ~ 0, .default = Other)) |>
  rename(Original = Other) |>
  select(time_value = Date, all_of(vocs)) |>
  mutate(
    total_count = rowSums(across(all_of(vocs))),
    across(all_of(vocs), ~ .x / total_count)
  ) |>
  select(-total_count)

# vocs "begin" when they first surpass 50% of sequences
voc_ranges <- us_variants |>
  pivot_longer(-time_value, names_to = "variant") |>
  summarise(start = min(time_value[value > 0.5]), .by = variant) |>
  mutate(
    end = pmin(
      lead(start) - 1,
      max(hhs$time_value, nchs$time_value),
      na.rm = TRUE
    )
  )

naive_hfr_offset <- 14L
geo_hfr_by_variant <- full_join(
  hhs,
  # move deaths backward to correspond with earlier hospitalizations
  nchs |> mutate(time_value = time_value - naive_hfr_offset),
  by = join_by(time_value, geo_value)
) |>
  left_join(voc_ranges, by = join_by(between(time_value, start, end))) |>
  select(-start, -end) |>
  summarise(
    hfr = sum(deaths, na.rm = TRUE) / sum(hosps, na.rm = TRUE),
    .by = c(geo_value, variant)
  )


daily_variants_all_geos <- bind_rows(
  us_variants |>
    filter(between(time_value, first_date, last_jhu)) |>
    complete(time_value = full_seq(c(first_date, last_jhu), period = 1)) |>
    mutate(
      across(all_of(vocs), ~ na.approx(.x, na.rm = FALSE)),
      geo_value = "us"
    ),
  # all states
  readRDS(here::here("Data", "Variants", "seq_prop_df.rds")) |>
    mutate(
      geo_value = tolower(State),
      Other = case_when(Date >= ymd("2022-01-01") ~ 0, .default = Other)
    ) |>
    rename(Original = Other) |>
    select(time_value = Date, geo_value, all_of(vocs)) |>
    mutate(
      total_prop = rowSums(across(all_of(vocs))),
      across(all_of(vocs), ~ .x / total_prop)
    ) |>
    select(-total_prop) |>
    arrange(geo_value, time_value) |>
    group_by(geo_value) |>
    complete(
      time_value = full_seq(c(time_value, first_date, last_jhu), period = 1)
    ) |>
    ungroup()
) |>
  group_by(geo_value) |>
  fill(all_of(vocs), .direction = "downup") |>
  ungroup() |>
  relocate(geo_value, time_value, all_of(vocs))

write_rds(
  daily_variants_all_geos,
  here::here("Data", "raw-data", "daily-variant-props.rds")
)

daily_hfrs <- daily_variants_all_geos |>
  pivot_longer(all_of(vocs), names_to = "variant", values_to = "prop") |>
  left_join(geo_hfr_by_variant, by = join_by(geo_value, variant)) |>
  mutate(hfr = hfr * prop) |>
  summarise(hfr = sum(hfr), .by = c(geo_value, time_value))

write_rds(
  daily_hfrs,
  here::here("Data", "raw-data", "daily-hfrs.rds")
)


# Compute Metadata for simulations --------------------------------------------

retro_meta_data <- summarise(
  full_retro,
  retro_lag_argmax = pmin(pmax(5, compute_max_acf_lag(hosps, deaths)), 35),
  .by = c(geo_value)
)

# For the real time case, we need to clean the reporting errors first
raw_finalized <- full_archive |> epix_as_of_current()
cleaned <- raw_finalized |>
  group_by(geo_value) |>
  mutate(
    deaths = replace_na(deaths, 0),
    hosps = replace_na(hosps, 0),
    ot = detect_outlr_rollmean(y = deaths, n = 42, detection_multiplier = 3)
  ) |>
  unpack(ot) |>
  mutate(
    deaths_no_outliers = case_when(!otlr ~ deaths, TRUE ~ repl),
    deaths_fix_reporting = fix_weekly_and_negatives(deaths_no_outliers, 1, 28),
    deaths_7dav = roll_meanr(deaths_fix_reporting, 7L),
    deaths_mu = roll_mean(deaths_fix_reporting, n = 28L, fill = NA)
  ) |>
  # remove NAs created by the above processing
  filter(between(time_value, min(time_value) + 15, max(time_value) - 15))

us_cleaned <- cleaned |>
  # outlier detection fails nationally (since outliers and nonoutliers are
  # aggregated together), so we recompute psuedo-national data by summing the
  # cleaned state-level data
  ungroup() |>
  filter(geo_value != "us") |>
  summarise(
    across(c(hosps, deaths, starts_with("deaths_")), sum),
    .by = "time_value"
  ) |>
  mutate(geo_value = "us")

cleaned <- bind_rows(cleaned |> filter(geo_value != "us"), us_cleaned)

write_rds(cleaned, here::here("Data", "raw-data", "cleaned-latest.rds"))

rt_meta_data <- cleaned |>
  summarise(
    dispersion = summary(
      glm(
        deaths_fix_reporting ~ deaths_mu + as.factor(wday(time_value)),
        family = quasipoisson()
      )
    )$dispersion,
    rt_lag_argmax = pmin(pmax(6, compute_max_acf_lag(hosps, deaths_7dav)), 35),
    rt_lag_1se = pmin(pmax(6, compute_1se_acf_lag(hosps, deaths_7dav)), 35)
  )

all_meta_data <- full_join(
  retro_meta_data,
  rt_meta_data,
  by = join_by(geo_value)
)

write_rds(
  all_meta_data,
  here::here("Data", "raw-data", "all-meta-data.rds")
)
