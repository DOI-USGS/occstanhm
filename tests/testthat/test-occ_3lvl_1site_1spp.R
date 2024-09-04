library(occstanhm)
library(tidyverse)

psi_sim <- 0.75
n_revisits <- 40
n_samples <- rpois(n_revisits, lambda = 10) + 1L
theta_sim <- 0.6
n_subsamples_in <- 10
p_sim <- 0.25

dat <-
  tibble(unit = "one",
         revisit = seq(1, n_revisits),
         z = rbinom(n = n_revisits, size = 1, prob = psi_sim),
         n_subsamples = n_subsamples_in) |>
  full_join(tibble(unit = "one",
                   revisit = rep(seq(1, n_revisits), times = n_samples),
                   theta =  theta_sim,
                   p = p_sim),
            by = c("unit", "revisit")) |>
  mutate(a = z * rbinom(n(), 1, theta),
         y = a * rbinom(n(), n_subsamples, p)) |>
  group_by(unit, revisit) |>
  mutate(a_obs = ifelse(a > 0, 1, 0),
         z_obs = ifelse(sum(y) > 0, 1, 0)) |>
  ungroup() |>
  rowid_to_column("index")

sample_summary <-
  dat |>
  group_by(unit, revisit) |>
  summarize(n_samples_per_revisit = n(),
            p_obs = mean(y / n_subsamples),
            theta_obs = mean(a_obs),
            z_obs = ifelse(sum(a_obs) > 0, 1, 0),
            unit_start = min(dplyr::cur_group_rows()),
            unit_stop = max(dplyr::cur_group_rows()),
            .groups = "drop")

stan_data <-
  list(y = dat |> pull(y),
       a_obs = dat |> pull(a_obs),
       k_subsamples = dat |> pull(n_subsamples),
       z_obs = sample_summary |> pull(z_obs),
       n_revisits = sample_summary |> nrow(),
       n_samples_per_revisit = sample_summary |> pull(n_samples_per_revisit),
       unit_start = sample_summary |> pull(unit_start),
       unit_stop = sample_summary |> pull(unit_stop),
       n_total_samples = dat |> nrow())

n_chains <- 1
n_parallel_chains <- 1
n_threads_per_chain <- 1
n_refresh <- 100
n_warmup <- 1
n_sample <- 1


suppressMessages({
  suppressWarnings({
    psi_theta_p_fit <-
      occ_3lvl_1site_1spp(stan_data,
                          n_chains = n_chains,
                          n_parallel_chains = n_parallel_chains,
                          n_threads_per_chain = n_threads_per_chain,
                          n_refresh = n_refresh,
                          n_warmup = n_warmup,
                          n_sample = n_sample)
  })
})

test_that("psi theta p model produces a  cmdstan fit", {
  expect_s3_class(psi_theta_p_fit, c("CmdStanMCMC", "CmdStanFit", "R6"))
})
