library(occstanhm)
library(tidyverse)
library(ggthemes)

## Set simulation details
# Number of units
n_units_sim <- 6 # max is 26 due to number of letters

# Number of species
n_species_sim <- 3

# Number of subunits (such as habitats within unit)
n_subunit <- 3

# number of markers
n_markers <- 2

# number of molecular samples
n_subsamples <- 6

# Poisson mean for number of stochastic samples per site
sample_lamba <- 10

## Simulate data
set.seed(12345)
# create species-unit data
dat_unit <-
  expand_grid(unit = letters[seq(1, n_units_sim)],
              species = paste("species", seq(1, n_species_sim))) |>
  mutate(psi = rbeta(n(), 6, 3))
dat_unit

# create subunit details
dat_subunit <-
  expand_grid(unit = letters[seq(1, n_units_sim)],
              subunit = paste("subunit", seq(1, n_subunit))) |>
  full_join(dat_unit, by = "unit", relationship = "many-to-many") |>
  arrange(unit, subunit, species) |>
  mutate(z = rbinom(n(), 1, psi),
         theta = rbeta(n(), 6, 3),
         unit_subunit = paste(unit, subunit, sep = "-"),
         subunit_species = paste(subunit, species, sep = "-"))

# create sample data
n_samples_per_subunit <-
  rpois(n = n_units_sim * n_subunit, lambda = sample_lamba) + 1L

dat_sample <-
  tibble(unit_subunit = rep(dat_subunit |>
                              distinct(unit_subunit) |>
                              pull(unit_subunit),
                            times = n_samples_per_subunit)) |>
  full_join(dat_subunit, by = "unit_subunit",
            relationship = "many-to-many") |>
  mutate(a = z * rbinom(n(), 1, theta))

# summarize data to see simulated theta values
theta_obs <-
  dat_sample |>
  group_by(unit, subunit, subunit_species, species, theta) |>
  summarize(n = n(), z = mean(z), theta_obs = mean(a),
            .groups = "drop") |>
  filter(z > 0) |>
  dplyr::select(unit, subunit_species, subunit, species, n, z, theta, theta_obs)

# Create/simulate full data set
# Arranged by unit, species, and then subunit
# Subunit arrange is tricky part (needs to be last)
dat <-
  expand_grid(unit = letters[seq(1, n_units_sim)],
              species = paste("species", seq(1, n_species_sim)),
              plain_marker = paste("M", seq(1, n_markers))) |>
  mutate(marker = paste(species, plain_marker, sep = "-"),
         unit_id = as.integer(factor(unit))) |>
#  select(-plain_marker) |>
  mutate(p = rbeta(n(), 3, 3)) |>
  full_join(dat_sample, by = c("unit", "species"),
            relationship = "many-to-many") |>
#  select(unit_id, unit, subunit, species, marker, z, a, p) |>
  mutate(subunit_species = paste(subunit, species, sep = "-"),
         k_subsamples = n_subsamples,
         y = a * rbinom(n(), k_subsamples, prob = p),
         a_obs = ifelse(y > 0, 1, 0)) |>
  group_by(unit_id, unit, subunit, species) |>
  mutate(z_obs = ifelse(sum(a_obs) > 0, 1, 0),
         subunit_id = as.integer(factor(subunit))) |>
  ungroup() |>
  arrange(unit, species, subunit) |>
  rowid_to_column("index")

# look at simulated or "observed" psi values
psi_obs <-
  dat |>
  group_by(unit_id, unit, species) |>
  filter(a > 0) |>
  summarize(p = mean(p), p_obs = mean(y / k_subsamples),
            .groups = "drop")

# Summarize simulated data for model inputs
# Need to make sure group by unit first
# This gives us  the start and stop for each species/unit combination
unit_species_summary <-
  dat |>
  group_by(unit_id, unit, species) |>
  summarize(psi_obs = mean(z_obs),
            unit_species_dat_start = min(index),
            unit_species_dat_stop = max(index),
            n_revisits_per_site = n(),
            .groups = "keep") |>
  mutate(logit_psi = qlogis(psi_obs))

# summarize each subunit or revisit to a unit
revisit_species_summary <-
  dat |>
  group_by(unit_id, unit, species, subunit) |>
  summarize(n_samples_per_grouping = n(),
            revisit_start = min(index),
            revisit_stop = max(index),
            z_obs = mean(z_obs),
            .groups = "drop") |>
  group_by(unit_id, unit, species) |>
  mutate(revisit_start_in_unit_species = (revisit_start -
                                            min(revisit_start) + 1L),
         revisit_stop_in_unit_species = revisit_stop - min(revisit_start) + 1L)

# summarize for slicing index
unit_species_revisit_summary <-
  revisit_species_summary |>
  rowid_to_column("index") |>
  group_by(unit_id, unit, species) |>
  summarize(unit_species_revisit_start = min(index),
            unit_species_revisit_stop = max(index),
            n_revisits_per_grouping = n(),
            psi_obs = mean(z_obs),
            .groups = "drop")

## Create parameters for model
# psi-level
x_psi <- model.matrix(~ species - 1, revisit_species_summary)

unit_species_summary_hyper <-
  unit_species_summary |>
  group_by(unit) |>
  summarize(n = n(),
            hyper_psi = mean(psi_obs),
            hyper_logit = mean(logit_psi))

x_psi_star <-
  model.matrix(~ 1, unit_species_summary_hyper)

beta_psi_star_prior <-
  matrix(0,
         nrow = x_psi_star |> ncol(),
         ncol = x_psi |> ncol())

# theta-level
w_theta <- model.matrix(~ subunit_species - 1, dat)
colnames(w_theta) <- gsub("subunit_species", "", colnames(w_theta))

w_theta_star <-
  model.matrix(~ 1, unit_species_summary_hyper)

alpha_theta_star_prior <-
  matrix(0,
         nrow = w_theta_star |> ncol(),
         ncol = w_theta |> ncol())

# p-level
v_p <- model.matrix(~ marker - 1, dat)
v_p_star <-
  model.matrix(~ 1, unit_species_summary_hyper)

delta_p_star_prior <-
  matrix(0,
         nrow = v_p_star |> ncol(),
         ncol = v_p |> ncol())

## Save data in list for Stan
stan_data <-
  # observation and latent observation data
  list(
    y = dat |> pull(y),
    k_subsamples = dat |> pull(k_subsamples),
    a_obs = dat |> pull(a_obs),
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

    # all different dimensions, starting with "n_"
    n_units = dat |> distinct(unit) |> nrow(),
    n_groupings = unit_species_summary |> nrow(),
    n_total_samples = dat |> nrow(),
    n_revisits_per_grouping = unit_species_revisit_summary |>
      pull(n_revisits_per_grouping),
    n_samples_per_grouping = revisit_species_summary |>
      pull(n_samples_per_grouping),
    n_total_groupings = revisit_species_summary |> nrow(),

    # psi-level inputs
    x_psi = x_psi,
    m_beta = x_psi |> ncol(),
    x_psi_star = x_psi_star,
    m_beta_star = x_psi_star |> ncol(),
    beta_psi_star_prior = beta_psi_star_prior,
    beta_psi_star_sd_prior = 1,
    eta_psi = 1,
    jj_psi = revisit_species_summary |> dplyr::pull(unit_id),

    # theta-level inputs
    w_theta = w_theta,
    m_alpha = w_theta |> ncol(),
    w_theta_star = w_theta_star,
    m_alpha_star = w_theta_star |> ncol(),
    alpha_theta_star_prior = alpha_theta_star_prior,
    alpha_theta_star_sd_prior = 1,
    eta_theta = 1,
    jj_theta = dat |> dplyr::pull(unit_id),

    # theta-level inputs
    v_p = v_p,
    m_delta = v_p |> ncol(),
    v_p_star = v_p_star,
    m_delta_star = v_p_star |> ncol(),
    delta_p_star_prior = delta_p_star_prior,
    delta_p_star_sd_prior = 1,
    eta_p = 1,
    jj_p = dat |> dplyr::pull(unit_id),

    # priors
    sigma_psi_prior = 1.5,
    sigma_psi_prior_sd = 0.05,
    sigma_theta_prior = 1.5,
    sigma_theta_prior_sd = 0.05,
    sigma_p_prior = 2.0,
    sigma_p_prior_sd = 0.05,

    # reduce_sum setting
    grainsize = 1
  )

## Set simulation details
# this includes within and among chain parallel settings
# (something that most Stan users may not be familiar with)
n_chains <- 1
n_parallel_chains <- 1
n_threads_per_chain <- 1
n_refresh <- 100

n_warmup <- 1
n_sample <- 1

suppressMessages({
  suppressWarnings({
    fit <- occstanhm_dev_3(stan_data,
                           n_chains = n_chains,
                           n_parallel_chains = n_parallel_chains,
                           n_threads_per_chain = n_threads_per_chain,
                           n_refresh = n_refresh,
                           n_warmup = n_warmup,
                           n_sample = n_sample)
  })
})

test_that("occstanhm_dev_3 produces a  cmdstan fit", {
  expect_s3_class(fit, c("CmdStanMCMC", "CmdStanFit", "R6"))
})
