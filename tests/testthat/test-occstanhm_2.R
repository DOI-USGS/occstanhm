library(occstanhm)
library(tidyverse)

stan_data <-
  system.file("extdata/occstanhm_2_stan_data.rds", package = "occstanhm") |>
  readRDS()

n_chains <- 1
n_parallel_chains <- 1
n_threads_per_chain <- 1
n_refresh <- 100
n_warmup <- 1
n_sample <- 1

suppressMessages({
  suppressWarnings({
    fit <-
      occstanhm_2(stan_data,
                  n_chains = n_chains,
                  n_parallel_chains = n_parallel_chains,
                  n_threads_per_chain = n_threads_per_chain,
                  n_refresh = n_refresh,
                  n_warmup = n_warmup,
                  n_sample = n_sample)
  })
})

test_that("occtan_2 produces a  cmdstan fit", {
  expect_s3_class(fit, c("CmdStanMCMC", "CmdStanFit", "R6"))
})
