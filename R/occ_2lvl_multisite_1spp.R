#' two-level single species occupancy model
#'
#' This is two-level occupancy model for a single species. The code is not
#' optimized and the file exists as a stepping stone to more complex models.
#' For example, a for loop would be replaced by a binomial logistic function
#' if the code were being written to be optimized.
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
#' n_units_sim <- 5 # max is 26 due to number of letters
#' psi_vec_in <- rbeta(n_units_sim, 3, 3)
#' theta_vec_in <- rbeta(n_units_sim, 3, 3)
#'
#' hist(psi_vec_in)
#' hist(theta_vec_in)
#'
#' n_units_in <- length(psi_vec_in)
#'
#' n_revisits_per_unit_in <- rpois(n_units_sim, 30)
#' hist(n_revisits_per_unit_in)
#'
#' n_samples_per_revisits <- rpois(sum(n_revisits_per_unit_in), 30) + 1L
#' hist(n_samples_per_revisits)
#'
#' ## unit and revisit data
#' dat_revisit <-
#'   tibble(
#'     unit = rep(letters[seq(1, n_units_in)], times = n_revisits_per_unit_in),
#'     psi_vec = rep(psi_vec_in, times = n_revisits_per_unit_in),
#'     theta_vec = rep(theta_vec_in, times = n_revisits_per_unit_in)) |>
#'   mutate(z =  rbinom(n(), size = 1, prob = psi_vec)) |>
#'   rowid_to_column("index") |>
#'   group_by(unit) |>
#'   mutate(revisit_id = dplyr::cur_group_rows() -
#'            min(dplyr::cur_group_rows()) + 1L) |>
#'   ungroup() |>
#'   select(unit, revisit_id, psi_vec, theta_vec, z)
#'
#' # sample
#' dat_sample <-
#'   tibble(unit =
#'            rep(dat_revisit |> pull(unit),
#'                times = n_samples_per_revisits),
#'          revisit_id =
#'            rep(dat_revisit |> pull(revisit_id),
#'                times = n_samples_per_revisits),
#'   ) |>
#'   rowid_to_column("index") |>
#'   group_by(unit) |>
#'   mutate(sample_id = dplyr::cur_group_rows() -
#'            min(dplyr::cur_group_rows()) + 1L) |>
#'   select(-index) |>
#'   ungroup()
#'
#' # merge together
#' dat <-
#'   dat_revisit |>
#'   full_join(dat_sample, by = c("unit", "revisit_id")) |>
#'   rowid_to_column("index") |>
#'   mutate(a = z * rbinom(n(), 1, theta_vec)) |>
#'   group_by(unit, revisit_id) |>
#'   mutate(z_obs = ifelse(sum(z) > 0, 1, 0)) |>
#'   ungroup()
#'
#' ## Need to make sure arrange by unit first
#' ## then within unit
#' unit_summary <-
#'   dat |>
#'   group_by(unit) |>
#'   summarize(psi_obs = mean(z_obs),
#'             unit_dat_start = min(index),
#'             unit_dat_stop = max(index),
#'             n_revisits_per_unit = n()) |>
#'   mutate(logit_psi = qlogis(psi_obs))
#'
#' unit_summary
#'
#' revisit_summary <-
#'   dat |>
#'   group_by(unit, revisit_id) |>
#'   summarize(n_samples_per_revisit = n(),
#'             revisit_start = min(index),
#'             revisit_stop = max(index),
#'             z_obs = mean(z_obs),
#'             .groups = "drop") |>
#'   group_by(unit) |>
#'   mutate(revisit_start_in_unit = revisit_start - min(revisit_start) + 1L,
#'          revisit_stop_in_unit = revisit_stop - min(revisit_start) + 1L)
#'
#' unit_revisit_summary <-
#'   revisit_summary |>
#'   rowid_to_column("index") |>
#'   group_by(unit) |>
#'   summarize(unit_revisit_start = min(index),
#'             unit_revisit_stop = max(index),
#'             n_revisits_per_unit = n())
#' stan_data <-
#'   list(
#'     a = dat |> pull(a),
#'     unit_dat_start = unit_summary |> pull(unit_dat_start),
#'     unit_dat_stop = unit_summary |> pull(unit_dat_stop),
#'     revisit_start_in_unit = revisit_summary |> pull(revisit_start_in_unit),
#'     revisit_stop_in_unit = revisit_summary |> pull(revisit_stop_in_unit),
#'     n_total_samples = dat |> nrow(),
#'     z_obs = revisit_summary |> pull(z_obs),
#'     n_samples_per_revisit = revisit_summary |> pull(n_samples_per_revisit),
#'     n_total_revisits = revisit_summary |> nrow(),
#'     unit_revisit_start = unit_revisit_summary |> pull(unit_revisit_start),
#'     unit_revisit_stop = unit_revisit_summary |> pull(unit_revisit_stop),
#'     n_revisits_per_unit = unit_revisit_summary |> pull(n_revisits_per_unit),
#'     n_units = unit_summary |> nrow(),
#'     grainsize = 1
#'   )
#'
#' n_chains <- 4
#' n_parallel_chains <- 4
#' n_threads_per_chain <- 4
#' n_refresh <- 100
#' n_warmup <- 1000
#' n_sample <- 1000
#'
#' fit <- occ_2lvl_multisite_1spp(stan_data,
#'                                n_chains = n_chains,
#'                                n_parallel_chains = n_parallel_chains,
#'                                n_threads_per_chain = n_threads_per_chain,
#'                                n_refresh = n_refresh,
#'                                n_warmup = n_warmup,
#'                                n_sample = n_sample)
#'
#' fit$summary() |>
#'   filter(grepl("logit_psi|logit_theta", variable)) |>
#'   select(variable, mean, median) |>
#'   mutate(mean = round(mean, 2),
#'          median = round(median, 2),
#'          prob = plogis(mean))
#'
#' unit_summary |> pull(psi_obs) |> qlogis() |> round(2)
#' unit_summary |> pull(psi_obs) |> round(2)
#'
#' theta_vec_in |> round(2)
#'}
occ_2lvl_multisite_1spp <-
  function(stan_data,
           n_chains = 2,
           n_parallel_chains = 2,
           n_threads_per_chain = 4,
           n_refresh = 100,
           n_warmup = 200,
           n_sample = 200,
           ...) {

    stan_file <-
      system.file("./stan_models/tutorial/occ_2lvl_multisite_1spp.stan",
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
