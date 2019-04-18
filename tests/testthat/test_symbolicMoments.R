#v1 Ist es konsistenter, wenn quote(0) anstatt 0 zurückgegeben wird?
#t1 multinomial testen

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
  # this example is based on [Lak+15], Appendix B
  
  cov <-  matrix(c(6, 2, 2, 4), nrow = 2)
  mean <- c(3, 4)
  order <- c(2, 1)
  
  # in this case we have beta = c(2, 1) and
  # A₁₁ = 1/2, A₁₂ = A₂₁ = 1, A₂₂ = 3
  
  # the following summands should appear in the result:
  # (A₁₁)₂*A₁₂ = (1/2 + 1)*(1/2)*1 = 0.75
  # 2*A₁₁*(A₁₂)₂ = 2*(1/2)*(1 + 1)*(1) = 2
  # (A₁₂)₃ = (1 + 2)*(1 + 1)*(1) = 6
  # (A₁₁)₂*A₂₂ = (1/2 + 1)*(1/2)*3 = 2.25
  # 2*A₁₁*A₁₂*A₂₂ = 2*(1/2)*1*3 = 3
  # (A₁₂)₂*A₂₂ = 3*(1 + 1)*(1) = 6
  # the sum is 20
  # multplied by beta[1]²*beta[2] = 4 leads to 80
  
  moment <- symbolicMoments(distribution = 'gamma', missingOrders = order, cov = cov, mean = mean)
  expect_identical(eval(moment[[1]]), 80)
  
  # test if 'gamma' works for more than one given order
  moments <- symbolicMoments(distribution = 'gamma', missingOrders = rbind(c(2,1), c(1,3)), cov = cov, mean = mean)
  expect_identical(lapply(moments, eval), list(80, 540))
})
