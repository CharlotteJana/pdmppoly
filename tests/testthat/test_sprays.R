#======== todo =================================================================

context("spray functions")

test_that("increase_arity works as expected", {
  library(spray)
  
  # add variables at the end:
  expect_true(lone(1,3) == increase_arity(lone(1,1), 2:3))
  
  # add variables at beginning and end:
  expect_true(lone(2,3) == increase_arity(lone(1,1), c(1,3)))
  
  # add two new variables side a side:
  expect_true(linear(c(1, 0, 0, 1)) == increase_arity(linear(c(1,1)), 2:3))
  
  # add a lot of new variables:
  expect_true(linear(c(0,0,1,0,1,0)) == increase_arity(linear(c(1,1)), c(1:2,4,6)))
  
  # the NULL polynomial:
  expect_true(0*lone(1,3) == increase_arity(0*lone(1,1), 2))
})