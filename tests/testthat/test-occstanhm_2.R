library(occstanhm)
library(tidyverse)

sim_data <-
  system.file("extdata/occstanhm_2_sim_data.rds", package = "occstanhm") |>
  readRDS()

revist_species_data <-
  sim_data$spp_obs |>
  ungroup() |>
  group_by(unit, species) |>
  mutate(min_idx = min(index) - 1) |>
  ungroup() |>
  group_by(unit, unit_id, species, species_id, revisit) |>
  summarize(z_obs = ifelse(sum(y) > 0, 1, 0),
            n_samples = n(),
            revisit_dat_start = min(index - min_idx),
            revisit_dat_stop = max(index - min_idx),
            .groups = "drop") |>
  rowid_to_column("index") |>
  ungroup()

unit_species_data <-
  sim_data$spp_obs |>
  group_by(unit, unit_id, species, species_id) |>
  summarize(psi_obs = mean(z_obs),
            n_samples = n(),
            unit_spp_dat_start = min(index),
            unit_spp_dat_stop = max(index),
            .groups = "drop")

revisit_species_at_unit_data  <-
  revist_species_data |>
  group_by(unit, species) |>
  summarize(n_revisits = n(),
            revisit_unit_start = min(index),
            revisit_unit_stop = max(index),
            .groups = "drop") |>
 ungroup()

x_psi <- model.matrix(~ species - 1, revist_species_data)
colnames(x_psi) <- gsub("species", "", colnames(x_psi))

x_psi_star <-
  sim_data$spp_obs |>
  distinct(unit) |>
  model.matrix(object = ~ 1)
colnames(x_psi_star) <- gsub("species", "", colnames(x_psi_star))
beta_psi_star_prior <-
  matrix(1.0,
         nrow = x_psi_star |> ncol(),
         ncol = x_psi |> ncol())
v_p <- model.matrix(~ species - 1, sim_data$spp_obs)
colnames(v_p) <- gsub("species", "", colnames(v_p))
## hyper parameters
v_p_star <-
  sim_data$spp_obs |>
  distinct(unit) |>
  model.matrix(object = ~ 1)
delta_p_star_prior <-
   matrix(1,
          nrow = v_p_star |> ncol(),
          ncol = v_p |> ncol())

stan_data <-
# first are observed data
list(y = sim_data$spp_obs |> pull(y),
     z_obs = revist_species_data |> pull(z_obs),

     # next are the dimensions for data (including vectors)
     n_unit_species = unit_species_data |> nrow(),
     n_total_revisits = revist_species_data |> nrow(),
     n_total_samples = sim_data$spp_obs |> nrow(),

     n_revisits = revisit_species_at_unit_data |> pull(n_revisits),
     n_samples = revist_species_data |> pull(n_samples),

     m_beta = x_psi |> ncol(),
     n_beta = sim_data$spp_obs |> distinct(unit) |> nrow(),

     m_beta_star = x_psi_star |> ncol(),
     n_beta_star = x_psi_star |> nrow(),

     m_delta = v_p |> ncol(),
     n_delta = sim_data$spp_obs |> distinct(unit_id) |> nrow(),

     m_delta_star = v_p_star |> ncol(),
     n_delta_star = v_p_star |> nrow(),

     # indexing
     unit_spp_dat_start = unit_species_data |> pull(unit_spp_dat_start),
     unit_spp_dat_stop = unit_species_data |> pull(unit_spp_dat_stop),

     revisit_dat_start = revist_species_data |> pull(revisit_dat_start),
     revisit_dat_stop = revist_species_data |> pull(revisit_dat_stop),

     revisit_unit_start =
        revisit_species_at_unit_data |> pull(revisit_unit_start),
     revisit_unit_stop =
        revisit_species_at_unit_data |> pull(revisit_unit_stop),

     # psi-level inputs
     x_psi = x_psi,
     x_psi_star = x_psi_star,

     # p-level inputs
     v_p = v_p,
     v_p_star = v_p_star,

     # psi-level priors
     beta_psi_star_prior = beta_psi_star_prior,
     beta_psi_star_sd_prior = 0.1,
     eta_psi = 1.0,
     jj_psi = revist_species_data |> dplyr::pull(unit_id),

     # p-level priors
     delta_p_star_prior = delta_p_star_prior,
     delta_p_star_sd_prior = 0.1,
     eta_p = 1.0,
     jj_p = sim_data$spp_obs |> dplyr::pull(unit_id),

     # priors for sigams
     sigma_psi_prior = 2.5,
     sigma_psi_prior_sd = 0.05,
     sigma_p_prior = 1.75,
     sigma_p_prior_sd = 0.05,

     # setting for reduce_sum
     grainsize = 1,

     ## Store data (not used in Stan)
     dat = sim_data$spp_obs,
     unit_species_data = unit_species_data
)

n_chains <- 1
n_parallel_chains <- 1
n_threads_per_chain <- 1
n_refresh <- 100
n_warmup <- 1
n_sample <- 1

suppressMessages({
  suppressWarnings({
    fit <-
      occstanhm_2(stan_data,
                  n_chains = n_chains,
                  n_parallel_chains = n_parallel_chains,
                  n_threads_per_chain = n_threads_per_chain,
                  n_refresh = n_refresh,
                  n_warmup = n_warmup,
                  n_sample = n_sample)
  })
})

test_that("occtan_2 produces a  cmdstan fit", {
  expect_s3_class(fit, c("CmdStanMCMC", "CmdStanFit", "R6"))
})
