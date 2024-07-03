library(occstanhm)
library(tidyverse)

n_units_sim <- 10 # max is 26 due to number of letters
psi_vec_in <- rbeta(n_units_sim, 3, 3)
n_units_in <- length(psi_vec_in)

n_revisits_per_site_in <- rpois(n_units_sim, 35)

dat <-
  tibble(
    unit = rep(letters[seq(1, n_units_in)], times = n_revisits_per_site_in),
    psi_vec = rep(psi_vec_in, times = n_revisits_per_site_in)
  ) |>
  mutate(z_obs =  rbinom(n(), size = 1, prob = psi_vec)) |>
  rowid_to_column("index") |>
  ungroup()

## Need to make sure arrange by unit first
## then within unit
unit_summary <-
  dat |>
  group_by(unit) |>
  summarize(psi_obs = mean(z_obs),
            unit_start = min(index),
            unit_stop = max(index),
            n_revisits_per_site = n()) |>
  mutate(logit_psi = qlogis(psi_obs))

stan_data <-
  list(
    z_obs = dat |> pull(z_obs),
    n_revisits = dat |> nrow(),
    n_revisits_per_unit = unit_summary |> pull(n_revisits_per_site),
    n_units = unit_summary |> nrow(),
    unit_start = unit_summary |> pull(unit_start),
    unit_stop = unit_summary |> pull(unit_stop),
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
    multipsi_fit <-
      occ_1lvl_multisite_1spp(stan_data = stan_data,
                              n_chains = n_chains,
                              n_threads_per_chain = n_threads_per_chain,
                              n_warmup = n_warmup,
                              n_sample = n_sample,
                              n_refresh = 0)
  })
})

test_that("occ_1lvl_multisite_1spp produces a  cmdstan fit", {
  expect_s3_class(multipsi_fit, c("CmdStanMCMC", "CmdStanFit", "R6"))
})
