#' multiple psi and multiple species model
#'
#' This model primarily exists to demonstrate how the occupancy models are
#' constructed. The model estimates multiple psi parameters. The model is a
#' building block for the multispecies model with correlation structure. This
#' model includes multiple locations and multiple species.
#'
#' @param stan_data Data formatted for Stan.
#' @return fit Model fit from cmdstanr.
#'
#' @export
#' @examples
#' # Increase the number of chains to 4
#' # Increase the number of parallel chains and threads based upon your CPUs
#' # For example, an 8 CPU machine might do 4 chains on 2 threads
#' # or 2 chains on 2 threads that are run back-to-back.
#' # Likewise, the number of warmup and samples to 100s for debugging and
#' # 1,000s or 10,000 for actual mode fits.
#' \dontrun{
#' library(tidyverse)
#'
#' n_units_sim <- 3 # max is 26 due to number of letters
#' n_species_sim <- 4
#' n_revisits_per_site_in <- rpois(n_units_sim, lambda = 35)
#'
#' n_revisits_per_site_in
#' n_revisits_per_site_in |> sum()
#'
#' unit_revisit_in <-
#'   tibble(unit = rep(letters[seq(1, n_units_sim)],
#'                     times = n_revisits_per_site_in)) |>
#'   group_by(unit) |>
#'   mutate(revisit = dplyr::cur_group_rows() -
#'           min(dplyr::cur_group_rows()) + 1L) |>
#'   ungroup()
#'
#' unit_revisit_in  |> print(n = Inf)
#'
#' dat <-
#'   expand_grid(
#'     unit = letters[seq(1, n_units_sim)],
#'     species = paste("species", seq(1, n_species_sim), sep = "-")) |>
#'   mutate(unit_id = as.integer(factor(unit)),
#'          psi = rbeta(n(), shape1 = 3, shape2 = 3)) |>
#'   full_join(unit_revisit_in, by = "unit",
#'             relationship = "many-to-many") |>
#'   mutate(z_obs =  rbinom(n(), size = 1, prob = psi)) |>
#'   arrange(unit, species) |>
#'   rowid_to_column("index") |>
#'   ungroup()
#'
#' ## Need to make sure arrange by unit first
#' ## then within unit
#' ## Also, make sure that dat is in full form (e.g., no missing rows)
#' unit_species_summary <-
#'   dat |>
#'   group_by(unit, unit_id, species) |>
#'   summarize(unit_species_start = min(index),
#'             unit_species_stop = max(index),
#'             n_revisits_per_site = n(),
#'             psi_obs = mean(z_obs),
#'             .groups = "drop") |>
#'   mutate(logit_psi = qlogis(psi_obs)) |>
#'   ungroup()
#'
#' unit_species_summary
#'
#' x_psi <- model.matrix(~ species - 1, unit_species_summary)
#' x_psi |> head()
#'
#' stan_data <-
#'   list(z_obs = dat |> dplyr::pull(z_obs),
#'        n_revisits = dat |> nrow(),
#'        n_revisits_per_site = unit_species_summary |>
#'                                 dplyr::pull(n_revisits_per_site),
#'        n_units_species = unit_species_summary |> nrow(),
#'        unit_start = unit_species_summary |> dplyr::pull(unit_species_start),
#'        unit_stop = unit_species_summary |> dplyr::pull(unit_species_stop),
#'        n_revisits_per_unit = unit_species_summary |>
#'                                 dplyr::pull(n_revisits_per_site),
#'        x_psi = x_psi,
#'        m_beta = x_psi |> ncol(),
#'        n_units = 3,
#'        jj_psi = unit_species_summary |> dplyr::pull(unit_id),
#'        grainsize = 1
#'   )
#'
#' n_chains <- 4
#' n_parallel_chains <- 4
#' n_threads_per_chain <- 4
#' n_refresh <- 100
#' n_warmup <- 2000
#' n_sample <- 2000
#'
#'
#' fit <- occ_1lvl_multispp(stan_data,
#'                          n_chains = n_chains,
#'                          n_parallel_chains = n_parallel_chains,
#'                          n_threads_per_chain = n_threads_per_chain,
#'                          n_refresh = n_refresh,
#'                          n_warmup = n_warmup,
#'                          n_sample = n_sample)
#'
#' fit$summary() |>
#'   filter(grepl("beta_psi\\[", variable)) |>
#'   print(n = Inf) |>
#'   select(variable, mean, median) |>
#'   mutate(obs = unit_species_summary |> arrange(species) |> pull(psi_obs) |>
#'            qlogis())
#'
#' unit_species_summary |> pull(psi_obs)
#'}
occ_1lvl_multispp <- function(stan_data,
                              n_chains = 2,
                              n_parallel_chains = 2,
                              n_threads_per_chain = 4,
                              n_refresh = 100,
                              n_warmup = 200,
                              n_sample = 200,
                              ...) {

  stan_file <-
    system.file("./stan_models/tutorial/occ_1lvl_multispp.stan",
                package = "occstanhm")

  mod <- cmdstanr::cmdstan_model(stan_file,
                                 cpp_options = list(stan_threads = TRUE))
  fit <- mod$sample(stan_data,
                    chains = n_chains,
                    parallel_chains = n_parallel_chains,
                    threads_per_chain = n_threads_per_chain,
                    refresh = n_refresh,
                    iter_warmup = n_warmup,
                    iter_sampling = n_sample,
                    ...)
  return(fit)
}
