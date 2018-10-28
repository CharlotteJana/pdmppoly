#======== todo =================================================================

library(spray)
context("helper functions")

test_that("increase_arity works as expected", {

  # add variables at the end:
  expect_true(lone(1,3) == increase_arity(lone(1,1), 2:3))
  
  # add variables at beginning and end:
  expect_true(lone(2,3) == increase_arity(lone(1,1), c(1,3)))
  
  # add two new variables side a side:
  expect_true(linear(c(1, 0, 0, 1)) == increase_arity(linear(c(1,1)), 2:3))
  
  # add a lot of new variables:
  expect_true(linear(c(0,0,1,0,1,0)) == increase_arity(linear(c(1,1)), c(1:2,4,6)))
  
  # the NULL polynomial:
  expect_true(0*lone(1,5) == increase_arity(0*lone(1,1), 2:5))
})

test_that("getIndex works as expected", {
  expect_identical(getIndex("c", letters), 3L)
  expect_equal(getIndex(2.0, c(1:5, 2)), c(2, 6))
  expect_equal(getIndex(4L, c(1:5, 4.0001)), 4)
  expect_equal(getIndex(1+1e-9, -1:2), 3)
  expect_error(getIndex(5, 1:3))
})

test_that("blowupSpray works as expected", {
  data("simplePoly")
  
  # polynomial without discrete variable
  s1 <- linear(1:0, 2)
  s2 <- linear(c(1,0,0,0), 2) 
  expect_true(blowupSpray(simplePoly, s1) == s2)
  
  # polynomial with only discrete variable
  s1 <- 3*lone(2, 2)
  s2 <- -3*lone(2, 4) + 3*lone(4, 4)
  expect_true(blowupSpray(simplePoly, s1) == s2)
  
  # polynomial with both variables
  s1 <- 5*product(2:3) + lone(1, 2)
  s2 <- -5*product(c(2,1,0,0)) + 5*product(c(2,0,0,1)) + lone(1, 4)
  expect_true(blowupSpray(simplePoly, s1) == s2)
  
  # with d not beeing the last variable in slot init
  init(simplePoly) <- c(d = 1, f = 2)
  s1 <- 5*product(3:2) + lone(2, 2)
  s2 <- -5*product(c(2,1,0,0)) + 5*product(c(2,0,0,1)) + lone(1, 4)
  expect_true(blowupSpray(simplePoly, s1) == s2)
})

test_that("blowupSpray works for model with 2 continous variables", {
  
  data("genePolyT")
  s1 <- 2*product(1:3) + lone(2, 3) + 3*lone(3,3)
  s2 <- 2*product(c(1,2,0,0,0,0))*(lone(3,6) + 2^3*lone(4,6) + 3^3*lone(5,6) + 4^3*lone(6,6)) + 
    lone(2,6) + 3*lone(3,6) + 6*lone(4,6) + 9*lone(5,6) + 12*lone(6,6)
  expect_true(blowupSpray(genePolyT, s1) == s2)
})