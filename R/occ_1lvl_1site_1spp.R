#' psi-only model
#'
#' This model primarily exists to demonstrate how the occupancy models are
#' constructed. The model estimates a site-level occupancy probability, psi,
#' on the logit scale with the output logit_psi.
#'
#' More practically, one would use a logistic regression. This model is included
#' to show the derivation of the model.
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
#' psi <- 0.25
#' n_revisits <- 400
#' z <- rbinom(n = n_revisits, size = 1, prob = psi)
#' n_chains <- 4
#' n_parallel_chains <- 4
#' n_threads_per_chain <- 4
#' n_refresh <- 100
#' n_warmup <- 2000
#' n_sample <- 2000
#'
#' stan_data <-
#'   list(z = z,
#'        n_revisits = n_revisits)
#'
#' psi_only_fit <-
#'   occ_1lvl_1site_1spp(stan_data = stan_data,
#'                       n_chains = n_chains,
#'                       n_threads_per_chain = n_threads_per_chain,
#'                       n_warmup = n_warmup,
#'                       n_sample = n_sample)
#' mean(z) |> qlogis() |> round(2)
#' mean(z)
#'
#'}
occ_1lvl_1site_1spp <- function(stan_data,
                                n_chains = 2,
                                n_parallel_chains = 2,
                                n_threads_per_chain = 4,
                                n_refresh = 100,
                                n_warmup = 200,
                                n_sample = 200,
                                ...) {

  stan_file <- system.file("./stan_models/tutorial/occ_1lvl_1site_1spp.stan",
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
