library(occstanhm)

n_units <- 50
n_unit_revisits <- 6
n_samples <- 10
n_subsamples <- 4
n_spp <- 4
spp_mu <- seq(-0.5, 1, length.out = n_spp)
sigma_psi <-
  sigma_theta <-
  sigma_p <-
  matrix(c(1.0,  0.0,  0.0,  0.0,
           0.0,  1.0,  1.0, -0.7,
           0.0,  1.0,  1.0, -0.7,
           0.0, -0.7, -0.7,  1.0) * 5,
         nrow = n_spp, ncol = n_spp)

example_out <-
  sim_long_form_occ_3(
    n_units = n_units,
    n_unit_revisit_mean = n_unit_revisits,
    k_samples_mean = n_samples,
    k_subsamples = n_subsamples,
    n_spp = n_spp,
    spp_mu = spp_mu,
    sigma_psi = sigma_psi,
    sigma_theta = sigma_psi,
    sigma_p = sigma_p
  )

test_that("length out outputs", {
  expect_equal(length(example_out), 4)
})

test_that("correlation inputs", {
  expect_equal(example_out$sigma_corr_psi$sigma,
               sigma_psi)
  expect_equal(example_out$sigma_corr_theta$sigma,
               sigma_theta)
  expect_equal(example_out$sigma_corr_p$sigma,
               sigma_p)
})

test_that("simulated matricies", {
  expect_equal(example_out$sigma_corr_psi$sigma_sim |> dim(),
               c(n_units, n_spp))
  expect_equal(example_out$sigma_corr_theta$sigma_sim |> dim(),
               c(n_units, n_spp))
  expect_equal(example_out$sigma_corr_p$sigma_sim |> dim(),
               c(n_units, n_spp))
})

test_that("simulated data", {
  expect_equal(example_out$sigma_corr_psi$unit_obs_spp |> dim(),
               c(n_units * n_spp, 4))
  expect_equal(example_out$sigma_corr_theta$unit_obs_spp |> dim(),
               c(n_units * n_spp, 4))
  expect_equal(example_out$sigma_corr_p$unit_obs_spp |> dim(),
               c(n_units * n_spp, 4))
})
