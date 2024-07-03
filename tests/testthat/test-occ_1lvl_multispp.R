library(occstanhm)
library(tidyverse)

n_units_sim <- 3 # max is 26 due to number of letters
n_species_sim <- 4
n_revisits_per_site_in <- rpois(n_units_sim, lambda = 35)

unit_revisit_in <-
  tibble(unit = rep(letters[seq(1, n_units_sim)],
                    times = n_revisits_per_site_in)) |>
  group_by(unit) |>
  mutate(revisit = dplyr::cur_group_rows() -
           min(dplyr::cur_group_rows()) + 1L) |>
  ungroup()

dat <-
  expand_grid(unit = letters[seq(1, n_units_sim)],
              species = paste("species", seq(1, n_species_sim), sep = "-")) |>
  mutate(unit_id = as.integer(factor(unit)),
         psi = rbeta(n(), shape1 = 3, shape2 = 3)) |>
  full_join(unit_revisit_in, by = "unit",
            relationship = "many-to-many") |>
  mutate(z_obs =  rbinom(n(), size = 1, prob = psi)) |>
  arrange(unit, species) |>
  rowid_to_column("index") |>
  ungroup()

## Need to make sure arrange by unit first
## then within unit
## Also, make sure that dat is in full form (e.g., no missing rows)
unit_species_summary <-
  dat |>
  group_by(unit, unit_id, species) |>
  summarize(unit_species_start = min(index),
            unit_species_stop = max(index),
            n_revisits_per_site = n(),
            psi_obs = mean(z_obs),
            .groups = "drop") |>
  mutate(logit_psi = qlogis(psi_obs)) |>
  ungroup()

x_psi <- model.matrix(~ species - 1, unit_species_summary)

stan_data <-
  list(
    z_obs = dat |> dplyr::pull(z_obs),
    n_revisits = dat |> nrow(),
    n_revisits_per_unit = unit_species_summary |>
      dplyr::pull(n_revisits_per_site),
    n_units_species = unit_species_summary |> nrow(),
    unit_start = unit_species_summary |> dplyr::pull(unit_species_start),
    unit_stop = unit_species_summary |> dplyr::pull(unit_species_stop),
    x_psi = x_psi,
    m_beta = x_psi |> ncol(),
    n_units = 3,
    jj_psi = unit_species_summary |> dplyr::pull(unit_id),
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
    multipsi_fit <- occ_1lvl_multispp(stan_data = stan_data,
                                      n_chains = n_chains,
                                      n_threads_per_chain = n_threads_per_chain,
                                      n_warmup = n_warmup,
                                      n_sample = n_sample,
                                      n_refresh = 0)
  })
})

test_that("occ_1lvl_multispp produces a  cmdstan fit", {
  expect_s3_class(multipsi_fit, c("CmdStanMCMC", "CmdStanFit", "R6"))
})
