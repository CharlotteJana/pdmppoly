context("moment closure")

test_that("closureCoefficient works for 1 continous variable",  {
  expect_identical(closureCoefficient(j = 2, m = 4, n = 5),
                   10*1 - 10*3  + 5*6)
})

test_that("closureCoefficient works for 2 continous variables", {
  expect_identical(closureCoefficient(j = 2:3, m = 3, n = 5:6), -400)
  expect_identical(closureCoefficient(j = 2:3, m = 4, n = 5:6), -400)
  expect_identical(closureCoefficient(j = 2:3, m = 5, n = 5:6), 0)
})

test_that("closureCoefficient throws errors", {
  expect_error(closureCoefficient(j = 1, m = 2, n = 5:6))
  expect_error(closuerCoefficient(j = 5, m = 2, n = 3))
  expect_error(closureCoefficient(j = 2, m = 5, n = 3))
})