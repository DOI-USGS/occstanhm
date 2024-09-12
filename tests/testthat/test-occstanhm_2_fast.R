library(occstanhm)
library(tidyverse)

n_chains <- 1
n_parallel_chains <- 1
n_threads_per_chain <- 1
n_refresh <- 100
n_warmup <- 1
n_sample <- 1

## @knitr simulate_data
n_units <- 50
n_unit_revisit_mean <- 10
n_unit_revisit_use <- 10
k_samples_mean <- 1
k_use <- 1
n_spp <- 3

## create covariance matrixes
a_sim_psi <- matrix(rbeta(n_spp^2, 0.01, 0.01) * 2 - 1, ncol = n_spp) 
sigma_use_psi <- t(a_sim_psi) %*% a_sim_psi
sigma_use_p <- matrix(0, nrow = n_spp, ncol = n_spp)
diag(sigma_use_p) <- 1

## create empty predictor matrices
x_psi_sim <-
  tibble(x_1 = rep(1, n_units), # intercept
         x_2 = rnorm(n_units, mean = 0.0, sd = 1.0),
         x_3 = rnorm(n_units, mean = 0.0, sd = 1.0))

v_p_sim <- NULL

psi_predictor_in <-
  rowSums(x_psi_sim)

sim_data <-
  sim_long_form_occ_2(n_units = n_units,
                      n_unit_revisit_mean = n_unit_revisit_mean,
                      n_spp = n_spp,
                      k_samples_mean = k_samples_mean,
                      spp_mu = rep(0.5, n_spp),
                      sigma_psi = sigma_use_psi,
                      sigma_p = sigma_use_p,
                      psi_predictor = psi_predictor_in,
                      p_predictor = NULL,
                      p_same = TRUE)

# known correlations among species
cor(sim_data$sigma_corr_psi$sigma_sim[,c(seq(1, n_spp, by = 1), n_spp)])

dat <-
  sim_data$spp_obs |>
  filter(sample <= k_use & revisit <= n_unit_revisit_use) |>
  mutate(species_id = as.integer(factor(species))) |>
  group_by(unit, unit_id, revisit, species, species_id) |>
  summarize(k_samples = n(),
            y = sum(y), .groups = "keep") |>
  mutate(z_obs = ifelse(sum(y) > 0, 1, 0)) |>
  ungroup()

clone_spp <-
  dat |>
  filter(species_id == max(species_id)) |>
  mutate(species_id = max(species_id) + 1L,
         species = paste0("Spp_", n_spp + 1))

dat <-
  dat |>
  bind_rows(clone_spp) |>
  arrange(unit_id, species_id) |>
  rowid_to_column("index")

revist_species_data <-
  dat |>
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
  ungroup() |>
  full_join(x_psi_sim |> rowid_to_column("unit"), by = "unit",
  )

revisit_species_at_unit_data  <-
  revist_species_data |>
  group_by(unit, species) |>
  summarize(n_revisits = n(),
            revisit_unit_start = min(index),
            revisit_unit_stop = max(index),
            .groups = "drop") |>
  ungroup()

unit_species_data <-
  dat |>
  group_by(unit, unit_id, species_id, species) |>
  summarize(psi_obs = mean(z_obs),
            n_samples = n(),
            unit_spp_dat_start = min(index),
            unit_spp_dat_stop = max(index),
            .groups = "drop")

x_psi <- model.matrix(~ species - 1 + x_1 + x_2 + x_3, revist_species_data)
colnames(x_psi) <- gsub("species", "", colnames(x_psi))

x_psi_star <-
  dat |>
  distinct(unit) |>
  model.matrix(object = ~ 1)

colnames(x_psi_star) <- gsub("species", "", colnames(x_psi_star))

## @knitr check_psi
x_psi_hyper <-
  dat |>
  group_by(unit, species) |>
  summarize(n = n(), .groups = "drop") |>
  mutate(one = 1)

## check math for beta coefficient
beta_psi_check <- matrix(1,
                         nrow = dat |> distinct(unit) |> nrow(),
                         ncol = x_psi |> ncol())

logit_psi_check <- matrix(1,
                          nrow = revist_species_data |> nrow(),
                          ncol = 1)

## @knitr beta_priors
## prior
beta_psi_star_prior <-
  matrix(0.0,
         nrow = x_psi_star |> ncol(),
         ncol = x_psi |> ncol())

## @knitr check_beta_priors
## check math for prior
x_psi_star_check <-
  matrix(1,
         nrow = x_psi_star |> nrow(),
         ncol = x_psi_star |> ncol())

## @knitr check_beta_star
beta_psi_star_check <-
  matrix(1,
         nrow = 1,
         ncol = dat |> distinct(species) |> nrow())

## @knitr p_coef
# p-level
v_p <- model.matrix(~ 1, dat)
colnames(v_p) <- gsub("species", "", colnames(v_p))

## @knitr p_coef_check
## check math for delta coefficient
delta_p_check <- matrix(1,
                        nrow = dat |> distinct(unit) |> nrow(),
                        ncol = v_p |> ncol())

logit_p_check <- matrix(1,
                        nrow = dat |> nrow(),
                        ncol = 1)

## @knitr save_as_stan
stan_data <-
  # first are observed data
  list(y = dat |> pull(y),
       k_samples = dat |> pull(k_samples),
       z_obs = revist_species_data |> pull(z_obs),
       
       # next are the dimensions for data (including vectors)
       n_unit_species = unit_species_data |> nrow(),
       n_total_revisits = revist_species_data |> nrow(),
       n_total_samples = dat |> nrow(),
       
       n_revisits = revisit_species_at_unit_data |> pull(n_revisits),
       n_samples = revist_species_data |> pull(n_samples),
       
       m_beta = x_psi |> ncol(),
       n_beta = dat |> distinct(unit) |> nrow(),
       
       m_beta_star = x_psi_star |> ncol(),
       n_beta_star = x_psi_star |> nrow(),
       
       m_delta = v_p |> ncol(),
       n_delta = 1, # dat |> distinct(unit_id) |> nrow(),
       
       # indexing
       unit_spp_dat_start = unit_species_data |> pull(unit_spp_dat_start),
       unit_spp_dat_stop = unit_species_data |> pull(unit_spp_dat_stop),
       
       revisit_dat_start = revist_species_data |> pull(revisit_dat_start),
       revisit_dat_stop = revist_species_data |> pull(revisit_dat_stop),
       
       revisit_unit_start = revisit_species_at_unit_data |> pull(revisit_unit_start),
       revisit_unit_stop = revisit_species_at_unit_data |> pull(revisit_unit_stop),
       
       # psi-level inputs
       x_psi = x_psi,
       x_psi_star = x_psi_star,
       
       # p-level inputs
       v_p = v_p,
       
       # psi-level priors
       beta_psi_star_prior = beta_psi_star_prior,
       beta_psi_star_sd_prior = 1,
       eta_psi = 0.1,
       jj_psi = revist_species_data |> dplyr::pull(unit_id),
       
       # p-level priors
       delta_p_star_sd_prior = 0.1,
       eta_p = 1.0,
       jj_p = rep(1, nrow(dat)), # dat |> dplyr::pull(unit_id),
       
       # priors for sigams
       sigma_psi_prior = 0.5,
       sigma_psi_prior_sd = 0.05,
       sigma_p_prior = 1.4,
       sigma_p_prior_sd = 0.05,
       
       # setting for reduce_sum
       grainsize = 1,
       
       ## Store data (not used in Stan)
       dat = dat,
       unit_species_data = unit_species_data
  )

## Look at correlation
psi_cor <-
  stan_data$dat |>
  dplyr::group_by(unit, species_id) |>
  dplyr::summarize(z_obs = mean(z_obs),
                   .groups = "drop") |>
  dplyr::select(unit, species_id, z_obs) |>
  pivot_wider(names_from = unit,
              values_from = z_obs) |>
  dplyr::select(-species_id) |>
  as.matrix() |>
  t() |>
  cor()

# get inputs in ball park of data
init_fun <- function(){
  list(logit_p = rnorm(n = stan_data$n_total_samples,
                       mean = qlogis(
                         (stan_data$y/stan_data$k_samples + .01) * 0.9),
                       sd = 0.05),
       delta_p = matrix(rnorm(1), nrow = 1, ncol = 1),
       logit_psi = rnorm(n = stan_data$n_total_samples,
                         mean = qlogis((stan_data$z + 0.01) * 0.9),
                         sd = 0.05),
       beta_psi = rnorm(n = n_units * (n_spp + 3),
                        mean = 0,
                        sd = 0.5) |>
         matrix(nrow = n_units, ncol = n_spp + 3, byrow = TRUE),
       l_omega_psi = (rbeta(n = choose(stan_data$m_beta, 2) - 1, 3, 2) * 2 - 1),
       tau_unif_psi = runif(stan_data$m_beta, min = 0, max = pi / 2),
       beta_psi_star = 
         matrix(rnorm(x_psi_star |> ncol() *  x_psi |> ncol(), 0, sd = 0.1),
                nrow = x_psi_star |> ncol(),
                ncol = x_psi |> ncol()),
       sigma_psi = abs(rnorm(1, 0, 1)),
       sigma_p = abs(rnorm(1, 0, 1)),
       z_psi =
         matrix(rnorm(stan_data$m_beta * stan_data$n_beta_star, 0, sd = 1),
                nrow = stan_data$m_beta,
                ncol = stan_data$n_beta_star),
       r_2_psi =
         runif(n = stan_data$m_beta - 1L, 0, 1)
  )}


suppressMessages({
  suppressWarnings({
    fit <-
      occstanhm_fast_2(stan_data,
                       n_chains = n_chains,
                       n_parallel_chains = n_parallel_chains,
                       n_threads_per_chain = n_threads_per_chain,
                       n_refresh = n_refresh,
                       n_warmup = n_warmup,
                       n_sample = n_sample)
  })
})

test_that("occtan_2_fast produces a  cmdstan fit", {
  expect_s3_class(fit, c("CmdStanMCMC", "CmdStanFit", "R6"))
})
