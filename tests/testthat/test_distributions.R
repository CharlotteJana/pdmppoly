context("distributions")

test_that("random.distribution returns a distribution", {
  d <- random.distribution(plot = FALSE, lower = 2, upper = 3)
  expect_true(is.list(d))
  expect_true(d[[1]]$spec %in% c("unif", "norm", "exp", "lnorm"))
})
