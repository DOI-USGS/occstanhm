#' Simulate a multivariate covariance/correlation matrix
#'
#' The `sigma` input allows you specify the default covariance. Leave
#' 1's along the diagonal to simulate according to a covariance matrix that is
#' equivalent to its correlation matrix. Alternatively,
#' you use random values, which is the default, by leaving this as `NULL`.
#'
#' The default is a three species site occupancy correlation matrix with 10
#' simulates units.
#'
#' The output include two formats, one as the raw matrix for use with functions
#' like `cor()` to examine the outputs and a second as tibble for use when
#' simulating data.
#'
#' @param n_units number of sampling units to simulate
#' @param n_spp number of species to simulate
#' @param spp_mu means for value
#' @param error_name name for MVN error
#' @param sigma_in MVN sigma, default is N
#' @export
#' @returns Returns a list with the following a list containing:
#'  \item{sigma_sim}{The simulated errors in matrix format.}
#'  \item{unit_obs_spp}{The simulated errors in tibble format.}
#'  \item{sigma}{The input value used for `sigma`}
#' @examples
#' # example code
#' example_out <-
#'   corr_sim(n_units = 10,
#'            n_spp = 3,
#'            spp_mu = c(0, 0, 0),
#'            error_name = "spp_error_psi",
#'            sigma_in = NULL)
#' cor(example_out$sigma_sim)

corr_sim <- function(
    n_units = 10,
    n_spp = 3,
    spp_mu = c(0, 0, 0),
    error_name = "spp_error_psi",
    sigma_in = NULL) {

  unit <- NULL

  unit_obs  <-
    tibble::tibble(unit = seq_len(n_units)) |>
    dplyr::mutate(unit_id = factor(unit))

  if (is.null(sigma_in)) {
    sigma_in <-
      matrix(stats::runif(n_spp^2) * 2 - 1, ncol = n_spp)
    sigma_in <-
      t(sigma_in) %*% sigma_in
  } else {
    if (all(dim(sigma_in) != c(n_spp, n_spp))) {
      stop("sigma_in must be NULL or an n_spp by n_spp matrix")
    }
  }

  sigma_sim <-
    matrix(mvtnorm::rmvnorm(n = n_units,
                            mean = spp_mu, sigma = sigma_in),
           ncol = n_spp, byrow = FALSE)
  colnames(sigma_sim) <- paste0("Spp_", seq_len(n_spp))
  sigma_sim <- as.data.frame(sigma_sim)

  unit_obs_spp <-
    dplyr::bind_cols(unit_obs, sigma_sim) |>
    tidyr::pivot_longer(cols = tidyr::starts_with("Spp_"),
                        names_to = "species",
                        values_to = error_name)

  return(
    list(
      sigma_sim = sigma_sim,
      unit_obs_spp = unit_obs_spp,
      sigma = sigma_in
    )
  )
}
