context("symbolic moments")

test_that("errors are thrown when wrong input", {
  
  cov <-  matrix(c(1, -1, 3, 2, 4, 5, 3, 5, 6), nrow = 3) # not symmetric
  expect_error(symbolicMoments(distribution = 'normal', cov = cov, missingOrders = 1:3))
  expect_identical(symbolicMoments(distribution = 'zero', cov = cov, missingOrders = 1:3), 
                   list(0)) # cov is not needed for 'zero' and therfore not validated
  
  cov[2] <- 2 # now it is symmetric
  expect_error(symbolicMoments(distribution = 'normal', cov = cov, missingOrders = 1:2)) # wrong dimensions of missingOrders
  expect_error(symbolicMoments(distribution = 'gamma', cov = cov, missingOrders = 1:3, mean = 1:2)) # wrong dimensions of mean
  
})

test_that("distribution = 'zero' works", {
  moment <- symbolicMoments(distribution = 'zero', missingOrders = matrix(1:6, nrow = 2))
  expect_equal(moment, list(0, 0))
})

test_that("distribution = 'normal' works", {
  
  cov <-  matrix(c("a", "b", "c", "d", 
                   "b", "e", "f", "h", 
                   "c", "f", "g", "i", 
                   "d", "h", "i", "j"), ncol = 4, byrow = TRUE)

  # order is odd -> moment should be zero
  moment <- symbolicMoments(distribution = 'normal', missingOrders = c(1,2,6,4), cov = cov)
  expect_equal(moment, list(0))
  
  # n = 4 and order is even (example: see Wikipedia on multivariate normal distribution)
  missingOrders <- matrix(c(4, 0, 0, 0,
                            3, 1, 0, 0,
                            2, 2, 0, 0,
                            2, 1, 1, 0,
                            1, 1, 1, 1), ncol = 4, byrow = TRUE)
  
  mom1 <- symbolicMoments(distribution = 'normal', 
                            missingOrders = missingOrders, cov = cov)
  mom2 <- list(quote(3 * "a"^2),
               quote(3 * ("a"^1 * "b"^1)),
               quote(2 * "b"^2 + 1 * ("a"^1 * "e"^1)),
               quote(2 * ("b"^1 * "c"^1) + 1 * ("a"^1 * "f"^1)),
               quote(1 * ("d"^1 * "f"^1) + 1 * ("c"^1 * "h"^1) + 1 * ("b"^1 * "i"^1)))
  expect_equal(mom1, mom2)
  
  # n = 1 and different orders
  mom1 <- symbolicMoments(distribution = 'normal', var = 4,
                          missingOrders = as.matrix(1:8, ncol = 1))
  mom1 <- sapply(mom1, eval)
  mom2 <- actuar::mnorm(1:8, mean = 0, sd = 2)
  expect_identical(mom1, mom2)
  
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
  skip("doesn't work")
  
  cov <-  matrix(c(6, 2, 2, 4), nrow = 2)
  mean <- c(3, 4)
  order <- c(2, 2)
  
  mom1 <- 4*log(cov[1,2] - mean[1]*mean[2]) +
          log(cov[1,1] - mean[1]^2) - 4*log(mean[1]) +
          log(cov[2,2] - mean[2]^2) - 4*log(mean[2])
  
  mom2 <- symbolicMoments(distribution = 'lognormal', missingOrders = order, cov = cov, mean = mean)
  expect_identical(mom1, eval(mom2[[1]]))
  
})

test_that("distribution = 'gamma' works for n = 2", {
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
