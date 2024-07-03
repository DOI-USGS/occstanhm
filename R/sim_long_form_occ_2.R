#' Simulate long-form 2-level occupancy data data
#'
#' This code simulates long-form (in contrast to wide, matrix format) data. The
#' function is used both to demonstrate and test the package. The function
#' simulates random numbers of revisits and samples to create more realistic
#' data, but also limit confidences that occurs when using fixed numbers for
#' these values.
#'
#' @param n_units Number of units to simulate
#' @param n_unit_revisit_mean Poisson average for revisits per unit
#' @param n_spp number of species
#' @param k_samples_mean Poisson average number of samples per revisit
#' @param spp_mu logit-scale average for p and psi for each species
#' @param sigma_psi variance/co-variance matrix for psi
#' @param sigma_p variance/co-variance matrix for p
#' @param psi_predictor prediction row (e.g., row sum of design matrix)
#' @param p_predictor prediction row (e.g., row sum of design matrix)
#' @param p_same Boolean if det prob p is the same across all species and units
#' @returns Returns a list with the following a list containing:
#'    \item{sigma_corr_psi}{The simulted error outputs used for psi}
#'    \item{sigma_corr_p}{The simulted error outputs used for p}
#'    \item{spp_obs}{A simulated dataset.}
#' @export
#' @examples
#' example_sim <-
#'   sim_long_form_occ_2()
sim_long_form_occ_2 <- function(
  n_units = 50,
  n_unit_revisit_mean = 6,
  n_spp = 4,
  k_samples_mean = 10,
  spp_mu = rep(0.5, n_spp),
  sigma_psi = matrix(c(1.0,  0.0,  0.0,  0.0,
                       0.0,  1.0,  1.0, -0.7,
                       0.0,  1.0,  1.0, -0.7,
                       0.0, -0.7, -0.7,  1.0),
                     nrow = n_spp, ncol = n_spp),
  sigma_p = matrix(c(1.0,  0.0,  0.0,  0.0,
                     0.0,  1.0,  1.0, -0.7,
                     0.0,  1.0,  1.0, -0.7,
                     0.0, -0.7, -0.7,  1.0),
                   nrow = n_spp, ncol = n_spp),
  psi_predictor = NULL,
  p_predictor = NULL,
  p_same = FALSE
) {

  spp_error_psi <- NULL
  spp_error_p <- NULL
  z_sim <- NULL
  species <- NULL
  y <- NULL
  unit <- NULL
  revisit <- NULL
  n_revisits_per_unit <- NULL
  n_samples <- NULL

  revisit_tibble <-
    tibble::tibble(unit = seq(1, n_units),
                   n_revisits_per_unit = rpois(n = n_units,
                                               n_unit_revisit_mean) + 1L) |>
    tidyr::uncount(n_revisits_per_unit) |>
    dplyr::group_by(unit) |>
    dplyr::mutate(revisit = seq(1, dplyr::n())) |>
    dplyr::ungroup()

  sample_tibble <-
    revisit_tibble |>
    dplyr::mutate(n_samples = rpois(n = dplyr::n(),
                                    lambda = k_samples_mean) + 1L) |>
    tidyr::uncount(n_samples) |>
    dplyr::group_by(unit, revisit) |>
    dplyr::mutate(sample = seq(1, dplyr::n())) |>
    dplyr::ungroup()

  sigma_corr_psi <-
    occstanhm::corr_sim(n_units = n_units,
                        n_spp = n_spp,
                        spp_mu = spp_mu,
                        sigma_in = sigma_psi,
                        error_name = "spp_error_psi")

  sigma_corr_p <-
    occstanhm::corr_sim(n_units = n_units, n_spp = n_spp,
                        spp_mu = spp_mu,
                        sigma_in = sigma_p,
                        error_name = "spp_error_p")
  if (p_same) {
    sigma_corr_p$unit_obs_spp <-
      sigma_corr_p$unit_obs_spp |>
      dplyr::mutate(spp_error_p = 0)
  }

  if (!is.null(psi_predictor)) {
    sigma_corr_psi$psi_predictor <- psi_predictor
    sigma_corr_psi$unit_obs_spp <-
      sigma_corr_psi$unit_obs_spp |>
      dplyr::mutate(spp_error_psi = spp_error_psi +
                      rep(psi_predictor, each = n_spp))
  }

  if (!is.null(p_predictor)) {
    sigma_corr_p$p_predictor <- p_predictor
    sigma_corr_p$unit_obs_spp <-
      sigma_corr_p$unit_obs_spp |>
      dplyr::mutate(spp_error_p = spp_error_p +
                      rep(p_predictor, each = n_spp))
  }

  spp_obs <-
    sigma_corr_psi$unit_obs_spp |>
    dplyr::full_join(sigma_corr_p$unit_obs_spp,
                     by = c("unit", "unit_id", "species")) |>
    dplyr::full_join(revisit_tibble, by = "unit",
                     relationship = "many-to-many") |>
    dplyr::ungroup() |>
    dplyr::mutate(z_sim = rbinom(n = dplyr::n(),
                                 size = 1,
                                 prob = plogis(spp_error_psi))) |>
    dplyr::full_join(sample_tibble, by = c("unit", "revisit"),
                     relationship = "many-to-many") |>
    dplyr::ungroup() |>
    dplyr::mutate(y = rbinom(n = dplyr::n(),
                             size = 1,
                             prob = plogis(spp_error_p)) *
                    z_sim) |>
    dplyr::group_by(unit, species, revisit) |>
    dplyr::mutate(z_obs = ifelse(sum(y) > 0, 1, 0)) |>
    dplyr::ungroup()

  return(list(sigma_corr_psi = sigma_corr_psi,
              sigma_corr_p = sigma_corr_p,
              spp_obs = spp_obs))

  return(sigma_corr_psi)
}
