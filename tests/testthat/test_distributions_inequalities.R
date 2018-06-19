#======== todo =================================================================
#t3 benötigte Pakete: distr <- wofür?

context("distributions - inequalities")


test_that("exists.distribution works as expected", {
  
  # beta distribution:
  expect_true(exists.distribution(0, 1, actuar::mbeta(1:4, 1, 2)))
  
  # uniform distribution:
  expect_true(exists.distribution(-1, 2, actuar::munif(1:4, min = -1, max = 2)))
  
  # bad moment values:
  expect_error(exists.distribution(-1, 2, -10))
  expect_false(exists.distribution(-1, 2, c(2, 0)))
  expect_false(exists.distribution(-1, 2, c(2, 5)))
  expect_false(exists.distribution(-4, -5, 1:5))
})


test_that("is.unimodal returns TRUE for simple unimodal distributions", {
  
  # uniform distribution:
  expect_true(is.unimodal(-10, 20, actuar::munif(1:4, min = -10, max = 20)))
  
  # beta distribution:
  expect_true(is.unimodal(0, 1, actuar::mbeta(1:4, 1, 2)))
})


test_that("is.unimodal returns FALSE for bimodal distributions", {
  
  # beta distribution:
  expect_false(is.unimodal(0, 1, actuar::mbeta(1:4, 0.5, 0.5)))
  
  # arcsine on [0,20]: 
  # see http://www.randomservices.org/random/special/Arcsine.html
  expect_false(is.unimodal(0, 20, c(20^1*1/2, 
                                    20^2*(1*3)/(2*4), 
                                    20^3*(1*3*5)/(2*4*6), 
                                    20^4*(1*3*5*7)/(2*4*6*8))))
})