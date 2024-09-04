library(occstanhm)
library(tidyverse)

psi <- 0.75
n_revisits <- 30
n_samples <- rpois(n_revisits, lambda = 10) + 1L
# +1 is to ensure at least 1 sample
theta_sim <- 0.6
dat <-
  tibble(unit = "one",
         revisit = seq(1, n_revisits),
         z = rbinom(n = n_revisits, size = 1, prob = psi)) |>
  full_join(tibble(unit = "one",
                   revisit = rep(seq(1, n_revisits), times = n_samples),
                   theta =  theta_sim),
            by = c("unit", "revisit")) |>
  mutate(a = z * rbinom(n(), 1, theta)) |>
  group_by(unit, revisit) |>
  mutate(z_obs = ifelse(sum(a) > 0, 1, 0)) |>
  ungroup() |>
  rowid_to_column("index")
sample_summary <-
  dat |>
  group_by(unit, revisit) |>
  summarize(n_samples_per_revisit = n(),
            theta_obs = mean(a),
            z_obs = ifelse(sum(a) > 0, 1, 0),
            unit_start = min(dplyr::cur_group_rows()),
            unit_stop = max(dplyr::cur_group_rows()),
            .groups = "drop")

stan_data <-
  list(a = dat |> pull(a),
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
    psi_theta_fit <-
      occ_2lvl_1site_1spp(stan_data = stan_data,
                          n_chains = n_chains,
                          n_threads_per_chain = n_threads_per_chain,
                          n_warmup = n_warmup,
                          n_sample = n_sample,
                          n_refresh = 0)
  })
})

test_that("occ_2lvl_1site_1spp model produces a  cmdstan fit", {
  expect_s3_class(psi_theta_fit, c("CmdStanMCMC", "CmdStanFit", "R6"))
})
