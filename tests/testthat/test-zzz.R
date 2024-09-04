test_that("cmdstan_version works correctly", {
  expect_no_error(occstanhm:::.onLoad(), message = NULL, class = NULL)
})
