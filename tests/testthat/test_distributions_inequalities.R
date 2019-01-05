context("is.unimodal")

test_that("is.unimodal returns `unimodal` for simple unimodal distributions", {
  
  # uniform distribution:
  expect_match(is.unimodal(-10, 20, actuar::munif(1:4, min = -10, max = 20)),
               "4-b-unimodal")
  
  # normal distribution with high values for a and b:
  expect_match(is.unimodal(-100, 100, actuar::mnorm(1:4)),
               "4-b-unimodal")
  
  # beta distribution:
  expect_match(is.4_b_unimodal(0, 1, actuar::mbeta(1:4, 1, 2)),
               "4-b-unimodal")
  expect_match(is.2_b_unimodal(0, 1, actuar::mbeta(1:2, 1, 2)),
               "2-b-unimodal")
})


test_that("is.unimodal returns `not unimodal` for bimodal distributions", {
  
  # beta distribution:
  expect_match(is.2_b_unimodal(0, 1, actuar::mbeta(1:2, 0.5, 0.5)), 
               "not unimodal")
  expect_match(is.4_b_unimodal(0, 1, actuar::mbeta(1:4, 0.5, 0.5)), 
               "not unimodal")
  
  # arcsine on [0,20]: 
  # see http://www.randomservices.org/random/special/Arcsine.html
  expect_match(is.unimodal(0, 20, c(20^1*1/2, 
                                    20^2*(1*3)/(2*4), 
                                    20^3*(1*3*5)/(2*4*6), 
                                    20^4*(1*3*5*7)/(2*4*6*8))),
               "not unimodal")
})

test_that("is.unimodal returns `not existant` in appropriate cases", {
  
  expect_error(is.unimodal(-1, 2, -10))
  expect_equal(is.unimodal(-1, 2, c(2, 5)), "not existant")
  expect_equal(is.unimodal(-4, -5, 1:5), "not existant")
})

test_that("is.unimodal works with vectorized input", {
  
  # moments given as matrix
  moments <- rbind(actuar::munif(1:2, min = -10, max = 10), actuar::mnorm(1:2))
  result <- is.unimodal(-10, 10, moments)
  expect_equal(result, c("2-b-unimodal", "2-b-unimodal"))
  
  # lower and upper given as vectors
  result <- is.unimodal(c(0, -10), c(1, -5), actuar::munif(1:4))
  expect_equal(result, c("4-b-unimodal", "not existant"))
  
  # moments as matrix, lower/upper as vectors with length = nrow(moments)
  moments <- rbind(actuar::mbeta(1:4, 0.5, 0.5), actuar::munif(1:4, min = -10, max = 10))
  result <- is.unimodal(c(0,-10), c(1,10), moments)
  expect_equal(result, c("not unimodal", "4-b-unimodal"))  
  
  # moments as matrix, lower/upper as vectors with length != nrow(moments)
  expect_error(is.unimodal(c(0,-10, 1), 10, moments))
})
