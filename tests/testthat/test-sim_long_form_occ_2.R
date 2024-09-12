library(occstanhm)

n_units <- 50
n_unit_revisit <- 6
n_spp <- 4
k_samples <- 10
spp_mu <- rep(0.5, n_spp)
sigma_psi <- matrix(c(1.0,  0.0,  0.0,  0.0,
                      0.0,  1.0,  1.0, -0.7,
                      0.0,  1.0,  1.0, -0.7,
                      0.0, -0.7, -0.7,  1.0),
                    nrow = n_spp, ncol = n_spp)
sigma_p <- matrix(c(1.0,  0.0,  0.0,  0.0,
                    0.0,  1.0,  1.0, -0.7,
                    0.0,  1.0,  1.0, -0.7,
                    0.0, -0.7, -0.7,  1.0),
                  nrow = n_spp, ncol = n_spp)

example_sim <-
  sim_long_form_occ_2(n_units = n_units,
                      n_unit_revisit = n_unit_revisit,
                      n_spp = n_spp,
                      k_samples = k_samples,
                      spp_mu = spp_mu,
                      sigma_psi = sigma_psi,
                      sigma_p = sigma_p)

test_that("sigam_corr_psi is correct", {
  expect_equal(length(example_sim$sigma_corr_psi),
               3)

  expect_equal(example_sim$sigma_corr_psi$sigma,
               sigma_psi)

  expect_equal(dim(example_sim$sigma_corr_psi$sigma_sim),
               c(n_units, n_spp))

  expect_equal(dim(example_sim$sigma_corr_psi$unit_obs_spp),
               c(n_units * n_spp, 4))
})

