#' Simulate long form data for a 3-level occupancy model
#'
#' Simulates data in a long form for a 3-level occupancy model
#
#' @param n_units Number of units to simulate
#' @param n_unit_revisit_mean Poisson average for revisits per unit
#' @param n_spp number of species
#' @param k_samples_mean Poisson average number of samples per revisit
#' @param k_subsamples number of subsamples per sample (molecular replicates)
#' @param spp_mu logit-scale average for p and psi for each species
#' @param sigma_psi variance/co-variance matrix for psi
#' @param sigma_theta variance/co-variance matrix for theta
#' @param sigma_p variance/co-variance matrix for p
#' @returns Returns a list with the following a list containing:
#'    \item{sigma_corr_psi}{The simulted error outputs used for psi}
#'    \item{sigma_corr_p}{The simulted error outputs used for p}
#'    \item{spp_obs}{A simulated dataset.}
#' @export
sim_long_form_occ_3 <- function(
  n_units = 50,
  n_unit_revisit_mean = 6,
  k_samples_mean = 10,
  k_subsamples = 4,
  n_spp = 4,
  spp_mu = seq(-0.5, 1, length.out = n_spp),
  sigma_psi = NULL,
  sigma_theta = NULL,
  sigma_p = NULL
) {
  species <- NULL
  revisit <- NULL
  z_sim <- NULL
  spp_error_psi <- NULL
  spp_error_theta <- NULL
  spp_error_p <- NULL
  a_sim <- NULL
  y <- NULL
  unit <- NULL
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
    corr_sim(n_units = n_units, n_spp = n_spp,
             spp_mu = spp_mu,
             sigma = sigma_psi, error_name = "spp_error_psi")
  sigma_corr_theta <-
    corr_sim(n_units = n_units, n_spp = n_spp,
             spp_mu = spp_mu,
             sigma = sigma_theta, error_name = "spp_error_theta")
  sigma_corr_p <-
    corr_sim(n_units = n_units, n_spp = n_spp,
             spp_mu = spp_mu,
             sigma = sigma_p, error_name = "spp_error_p")

  spp_obs <-
    sigma_corr_psi$unit_obs_spp |>
    dplyr::full_join(sigma_corr_theta$unit_obs_spp,
                     by = c("unit", "unit_id", "species")) |>
    dplyr::full_join(sigma_corr_p$unit_obs_spp,
                     by = c("unit", "unit_id", "species")) |>
    dplyr::full_join(revisit_tibble, by = "unit",
                     relationship = "many-to-many") |>
    dplyr::ungroup() |>
    dplyr::mutate(z_sim = rbinom(n = dplyr::n(),
                                 size = 1,
                                 prob = plogis(spp_error_psi)),
                  k_subsample = k_subsamples) |>
    dplyr::full_join(sample_tibble, by = c("unit", "revisit"),
                     relationship = "many-to-many") |>
    dplyr::ungroup() |>
    dplyr::mutate(
      a_sim = rbinom(n = dplyr::n(),
                     size = 1,
                     prob = plogis(spp_error_theta)) * z_sim,
      y = rbinom(n = dplyr::n(),
                 size = k_subsamples,
                 prob = plogis(spp_error_p)) * a_sim,
      a_obs = ifelse(y > 0, 1, 0)
    ) |>
    dplyr::group_by(unit, species, revisit) |>
    dplyr::mutate(z_obs = ifelse(sum(y) > 0, 1, 0)) |>
    dplyr::ungroup()

  return(list(sigma_corr_psi = sigma_corr_psi,
              sigma_corr_theta = sigma_corr_theta,
              sigma_corr_p = sigma_corr_p,
              spp_obs = spp_obs))
}
