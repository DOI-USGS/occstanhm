#' psi, theta, and p model
#'
#' This model primarily exists to demonstrate how the occupancy models are
#' constructed. The model estimates a site-level occupancy probability, psi, and
#' p on the logit scale with the.
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
#' psi_sim <- 0.75
#' n_revisits <- 40
#' n_samples <- rpois(n_revisits, lambda = 10) + 1L
#' # +1 is to ensure at least 1
#' n_samples
#' theta_sim <- 0.6
#' theta_sim
#' n_subsamples_in <- 10
#' p_sim <- 0.25
#'
#' dat <-
#'   tibble(unit = "one",
#'          revisit = seq(1, n_revisits),
#'          z = rbinom(n = n_revisits, size = 1, prob = psi_sim),
#'          n_subsamples = n_subsamples_in) |>
#'   full_join(tibble(unit = "one",
#'                    revisit = rep(seq(1, n_revisits), times = n_samples),
#'                    theta =  theta_sim,
#'                    p = p_sim),
#'             by = c("unit", "revisit")) |>
#'   mutate(a = z * rbinom(n(), 1, theta),
#'          y = a * rbinom(n(), n_subsamples, p)) |>
#'   group_by(unit, revisit) |>
#'   mutate(a_obs = ifelse(a > 0, 1, 0),
#'          z_obs = ifelse(sum(y) > 0, 1, 0)) |>
#'   ungroup() |>
#'   rowid_to_column("index")
#'
#' dat
#'
#' dat |>
#'   filter(z != z_obs)
#'
#' dat |>
#'   filter(a > 0)
#'
#' dat |>
#'   filter(a != a_obs)
#'
#' sample_summary <-
#'   dat |>
#'   group_by(unit, revisit) |>
#'   summarize(n_samples_per_revisit = n(),
#'             p_obs = mean(y/n_subsamples),
#'             theta_obs = mean(a_obs),
#'             z_obs = ifelse(sum(a_obs) > 0, 1, 0),
#'             unit_start = min(dplyr::cur_group_rows()),
#'             unit_stop = max(dplyr::cur_group_rows()),
#'             .groups = "drop")
#'
#' sample_summary
#'
#' stan_data <-
#'   list(y = dat |> pull(y),
#'        a_obs = dat |> pull(a_obs),
#'        k_subsamples = dat |> pull(n_subsamples),
#'        z_obs = sample_summary |> pull(z_obs),
#'        n_revisits = sample_summary |> nrow(),
#'        n_samples_per_revisit = sample_summary |> pull(n_samples_per_revisit),
#'        unit_start = sample_summary |> pull(unit_start),
#'        unit_stop = sample_summary |> pull(unit_stop),
#'        n_total_samples = dat |> nrow())
#'
#' n_chains <- 4
#' n_parallel_chains <- 4
#' n_threads_per_chain <- 4
#' n_refresh <- 100
#' n_warmup <- 2000
#' n_sample <- 2000
#'
#' fit <- occ_3lvl_1site_1spp(stan_data,
#'                            n_chains = n_chains,
#'                            n_parallel_chains = n_parallel_chains,
#'                            n_threads_per_chain = n_threads_per_chain,
#'                            n_refresh = n_refresh,
#'                            n_warmup = n_warmup,
#'                            n_sample = n_sample)
#'
#' fit_summary <-
#'   fit$summary()
#'
#' fit_summary |>
#'   filter(variable %in% c("logit_psi", "logit_theta", "logit_p")) |>
#'   select(variable, mean, median) |>
#'   mutate(p_mean = plogis(median))
#'
#' psi_sim
#' theta_sim
#' p_sim

#' sample_summary |>
#'   summarize(obs_psi = mean(z_obs)) |>
#'   mutate(logit_psi = qlogis(obs_psi))
#'}
occ_3lvl_1site_1spp <- function(stan_data,
                                n_chains = 2,
                                n_parallel_chains = 2,
                                n_threads_per_chain = 4,
                                n_refresh = 100,
                                n_warmup = 200,
                                n_sample = 200,
                                ...) {

  stan_file <- system.file("./stan_models/tutorial/occ_3lvl_1site_1spp.stan",
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
