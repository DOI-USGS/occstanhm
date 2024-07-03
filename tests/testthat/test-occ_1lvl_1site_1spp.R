library(occstanhm)

psi <- 0.25
n_revisits <- 400
z <- rbinom(n = n_revisits, size = 1, prob = psi)
n_chains <- 1
n_parallel_chains <- 1
n_threads_per_chain <- 1
n_refresh <- 100
n_warmup <- 1
n_sample <- 1

stan_data <-
  list(z = z,
       n_revisits = n_revisits)

suppressMessages({
  suppressWarnings({
    psi_fit <- occ_1lvl_1site_1spp(stan_data = stan_data,
                                   n_chains = n_chains,
                                   n_threads_per_chain = n_threads_per_chain,
                                   n_warmup = n_warmup,
                                   n_sample = n_sample,
                                   n_refresh = 0)
  })
})

test_that("occ_1lvl_1site_1spp produces a  cmdstan fit", {
  expect_s3_class(psi_fit, c("CmdStanMCMC", "CmdStanFit", "R6"))
})
