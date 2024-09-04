library(occstanhm)
library(tidyverse)

n_units_sim <- 2 # max is 26 due to number of letters
n_species_sim <- 3

n_revisits_per_site_in <- rpois(n_units_sim, 30)
n_samples_per_revisits <- rpois(sum(n_revisits_per_site_in), 20) + 1L

## species unit
dat_species <-
  expand_grid(unit = rep(letters[seq(1, n_units_sim)]),
              species = paste("species", seq(1, n_species_sim))) |>
  mutate(psi_vec = rbeta(n(), 3, 3),
         theta_vec = rbeta(n(), 3, 3),
         p_vec = rbeta(n(), 3, 3),
         k_subsamples = 6)

## unit and revisit data
dat_revisit <-
  tibble(unit = rep(letters[seq(1, n_units_sim)],
                    times = n_revisits_per_site_in)) |>
  rowid_to_column("index") |>
  group_by(unit) |>
  mutate(revisit_id = dplyr::cur_group_rows() -
           min(dplyr::cur_group_rows()) + 1L) |>
  ungroup() |>
  select(unit, revisit_id)

# sample
dat_sample <-
  tibble(
    unit = rep(dat_revisit |> pull(unit),
               times = n_samples_per_revisits),
    revisit_id = rep(dat_revisit |> pull(revisit_id),
                     times = n_samples_per_revisits),
  ) |>
  rowid_to_column("index") |>
  group_by(unit) |>
  mutate(sample_id = dplyr::cur_group_rows() -
           min(dplyr::cur_group_rows()) + 1L) |>
  select(-index) |>
  ungroup()

dat <-
  dat_species |>
  full_join(dat_revisit, by = "unit",
            relationship = "many-to-many") |>
  mutate(z = rbinom(n(), 1, psi_vec)) |>
  full_join(dat_sample, by = c("unit", "revisit_id"),
            relationship = "many-to-many") |>
  arrange(unit, species) |>
  rowid_to_column("index") |>
  mutate(a = z * rbinom(n(), 1, theta_vec),
         y = a * rbinom(n(), k_subsamples, p_vec),
         a_obs = ifelse(y > 0, 1, 0)) |>
  ungroup() |>
  group_by(unit, revisit_id, species) |>
  mutate(z_obs = ifelse(sum(a_obs) > 0, 1, 0)) |>
  ungroup()

unit_species_summary <-
  dat |>
  group_by(unit, species) |>
  summarize(psi_obs = mean(z_obs),
            unit_species_dat_start = min(index),
            unit_species_dat_stop = max(index),
            n_revisits_per_site = n(),
            .groups = "keep") |>
  mutate(logit_psi = qlogis(psi_obs))

revisit_species_summary <-
  dat |>
  group_by(unit, species, revisit_id) |>
  summarize(n_samples_per_revisit = n(),
            revisit_start = min(index),
            revisit_stop = max(index),
            z_obs = mean(z_obs),
            .groups = "drop") |>
  group_by(unit, species) |>
  mutate(revisit_start_in_unit_species = (revisit_start -
                                            min(revisit_start) + 1L),
         revisit_stop_in_unit_species = revisit_stop - min(revisit_start) + 1L)

unit_species_revisit_summary <-
  revisit_species_summary |>
  rowid_to_column("index") |>
  group_by(unit, species) |>
  summarize(unit_species_revisit_start = min(index),
            unit_species_revisit_stop = max(index),
            n_revisits_per_unit_species = n(),
            psi_obs = mean(z_obs),
            .groups = "drop")

stan_data <-
  # observation and latent observation data
  list(
    y = dat |> pull(y),
    k_subsamples = dat |> pull(k_subsamples),
    a_obs = dat |> pull(a),
    z_obs = revisit_species_summary |> pull(z_obs),
    # indexing
    unit_species_dat_start = unit_species_summary |>
      pull(unit_species_dat_start),
    unit_species_dat_stop = unit_species_summary |>
      pull(unit_species_dat_stop),
    revisit_start_in_unit = revisit_species_summary |>
      pull(revisit_start_in_unit_species),
    revisit_stop_in_unit = revisit_species_summary |>
      pull(revisit_stop_in_unit_species),
    unit_species_revisit_start = unit_species_revisit_summary |>
      pull(unit_species_revisit_start),
    unit_species_revisit_stop = unit_species_revisit_summary |>
      pull(unit_species_revisit_stop),
    # all different dimensions, such as "n"
    n_units_species = unit_species_summary |> nrow(),
    n_total_samples = dat |> nrow(),
    n_revisits_per_unit_species = unit_species_revisit_summary |>
      pull(n_revisits_per_unit_species),
    n_samples_per_revisit_species = revisit_species_summary |>
      pull(n_samples_per_revisit),
    n_total_revisits_species = revisit_species_summary |> nrow(),
    # reduce_sum setting
    grainsize = 1
  )

n_chains <- 1
n_parallel_chains <- 1
n_threads_per_chain <- 1
n_refresh <- 100
n_warmup <- 1
n_sample <- 1

suppressMessages({
  suppressWarnings({
    fit <- occ_3lvl_multispp(stan_data,
                             n_chains = n_chains,
                             n_parallel_chains = n_parallel_chains,
                             n_threads_per_chain = n_threads_per_chain,
                             n_refresh = n_refresh,
                             n_warmup = n_warmup,
                             n_sample = n_sample)
  })
})

test_that("occ_3lvl_multispp produces a  cmdstan fit", {
  expect_s3_class(fit, c("CmdStanMCMC", "CmdStanFit", "R6"))
})
