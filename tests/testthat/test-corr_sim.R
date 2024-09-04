library(occstanhm)

set.seed(1234)
n_spp <- 3
n_units <- 1000
sigma_in <-
  matrix(c(1,   0.1,   1,
           0.1,   1,   0,
           1,   0,   1),
         ncol = n_spp)

example_out <-
  suppressWarnings({
    corr_sim(n_units = n_units,
             n_spp = n_spp,
             spp_mu = c(0, 0, 0),
             error_name = "spp_error_psi",
             sigma = sigma_in)
  })

test_that("check input matrix", {
  expect_equal(example_out$sigma, sigma_in)
})

test_that("check simulated correlation matrix", {
  expect_equal(as.numeric(cor(example_out$sigma_sim)),
               as.numeric(sigma_in),
               tolerance = 0.05)
})

test_that("sigma set correctly", {
  expect_equal(example_out$sigma, sigma_in)
})

test_that("unit_obs_spp is correct dim", {
  expect_equal(example_out$unit_obs_spp |> dim(),
               c(n_spp * n_units, 4))
})

test_that("no error occurs if is not n_spp by n_spp", {
  expect_no_error(corr_sim(), message = NULL, class = NULL)
})

test_that("error occurs if sigma not n_spp by n_spp",
          {
            expect_error(corr_sim(n_units = 1000,
                                  n_spp = 200,
                                  spp_mu = c(0, 0, 0),
                                  error_name = "spp_error_psi",
                                  sigma = sigma_in),
                         "sigma_in must be NULL or an n_spp by n_spp matrix")
          })
