#v1 Ist es konsistenter, wenn quote(0) anstatt 0 zurückgegeben wird?

context("symbolic moments")

test_that("covariance matrix is validated", {
  cov <-  matrix(c(1, 2, 3, 2, 4, 5, 3, 5, 6), nrow = 3)
})

test_that("mean is validated", {
})

test_that("distribution = 'zero' works", {
  moment <- symbolicMoments(distribution = 'zero', missingOrders = matrix(1:6, nrow = 2))
  expect_equal(moment, list(0, 0))
})

test_that("distribution = 'normal' works", {
  skip("not implemented yet")
  
  cov <-  matrix(c(1, 2, 3, 2, 4, 5, 3, 5, 6), nrow = 3)
  
  # order is odd -> moment should be zero
  moment <- symbolicMoments(distribution = 'normal', missingOrders = c(1,2,4), cov = cov)
  expect_equal(moment, list(0))
  
  # order is even
  moment <- symbolicMoments(distribution = 'normal', missingOrders = 1:3, cov = cov)
  moment <- lapply(moment, eval)
})

test_that("distribution = 'lognormal' works for n = 1", {

  order <- 5
  meanlog <- 1
  sdlog <- 2
  mom1 <- actuar::mlnorm(order, meanlog = meanlog, sdlog = sdlog)
  
  mean <- exp(meanlog+0.5*sdlog^2)
  var <- exp(2*meanlog+sdlog^2)*(exp(sdlog^2)-1)
  symbolicMoments(distribution = 'lognormal', missingOrders = order, var = var, mean = mean)
  mom2 <- symbolicMoments(distribution = 'lognormal', missingOrders = order, var = var, mean = mean)
  mom2 <- sapply(mom2, eval)
  
  expect_equal(mom1, mom2)
})

test_that("distribution = 'lognormal' works for n = 2", {
  skip("not implemented yet")
})

test_that("distribution = 'gamma' works", {
  
  cov <-  matrix(c(1, 2, 3, 2, 4, 5, 3, 5, 6), nrow = 3)
  mean <- c("m1", "m2", "m3")
  moment <- symbolicMoments(distribution = 'gamma', missingOrders = 1:3, cov = cov, mean = mean)

})
