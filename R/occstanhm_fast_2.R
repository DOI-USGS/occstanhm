#' Two-level, multispecies hierarchical occupancy model
#'
#' occstanhm_2 is a two-level, hierarchical occupancy model.
#' occstanhm_fast_2 is a stripped down version of the two-level model designed
#' to be comparable to Tolber et al. 2019 (https://doi.org/10.1002/ecy.2754).
#'
#' The model does not included hyper-parameters and a correlation matrix for the
#' p-level fo the model. The model also allows for a single detection
#' probability to be estimated across all species.
#'
#' For an example, see the comarpison vignette.
#'
#' @param stan_data Data formatted for Stan.
#' @return fit Model fit from cmdstanr.
#'
#' @export
#' @examples
#'
#' \dontrun{
#'
#' }

occstanhm_fast_2 <- function(
    stan_data,
    n_chains = 2,
    n_parallel_chains = 2,
    n_threads_per_chain = 4,
    n_refresh = 100,
    n_warmup = 200,
    n_sample = 200,
    ...) {

  stan_file <- system.file("./stan_models/occstanhm_fast_2.stan",
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
