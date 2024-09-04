#' multiple psi model
#'
#' This model primarily exists to demonstrate how the occupancy models are
#' constructed. The model estimates multiple psi parameters. The model is a
#' building block for the multispecies model with correlation structure.
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
#' n_units_sim <- 10 # max is 26 due to number of letters
#' psi_vec_in <- rbeta(n_units_sim, 3, 3)
#' hist(psi_vec_in)
#' n_units_in <- length(psi_vec_in)
#'
#' n_revisits_per_unit_in <- rpois(n_units_sim, 35)
#'
#' dat <-
#'   tibble(
#'     unit = rep(letters[seq(1, n_units_in)], times = n_revisits_per_unit_in),
#'     psi_vec = rep(psi_vec_in, times = n_revisits_per_unit_in)) |>
#'   mutate(z_obs =  rbinom(n(), size = 1, prob = psi_vec)) |>
#'   rowid_to_column("index") |>
#'   ungroup()
#'
#' ## Need to make sure arrange by unit first
#' ## then within unit
#' unit_summary <-
#'   dat |>
#'   group_by(unit) |>
#'   summarize(psi_obs = mean(z_obs),
#'             unit_start = min(index),
#'             unit_stop = max(index),
#'             n_revisits_per_unit = n()) |>
#'   mutate(logit_psi = qlogis(psi_obs))
#'
#' unit_summary
#'
#' stan_data <-
#'   list(z_obs = dat |> pull(z_obs),
#'        n_revisits = dat |> nrow(),
#'        n_revisits_per_unit = unit_summary |> pull(n_revisits_per_unit),
#'        n_units = unit_summary |> nrow(),
#'        unit_start = unit_summary |> pull(unit_start),
#'        unit_stop = unit_summary |> pull(unit_stop),
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
#' fit <- occ_1lvl_multisite_1spp(stan_data = stan_data,
#'                                n_chains = n_chains,
#'                                n_parallel_chains = n_parallel_chains,
#'                                n_threads_per_chain = n_threads_per_chain,
#'                                n_refresh = n_refresh,
#'                                n_warmup = n_warmup,
#'                                n_sample = n_sample)
#'
#' fit
#' unit_summary |> pull(psi_obs) |> qlogis() |> round(2)
#' unit_summary |> pull(psi_obs)
#'
#'}
occ_1lvl_multisite_1spp <- function(stan_data,
                                    n_chains = 2,
                                    n_parallel_chains = 2,
                                    n_threads_per_chain = 4,
                                    n_refresh = 100,
                                    n_warmup = 200,
                                    n_sample = 200,
                                    ...) {

  stan_file <-
    system.file("./stan_models/tutorial/occ_1lvl_multisite_1spp.stan",
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
