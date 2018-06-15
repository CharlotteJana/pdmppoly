#======== todo =================================================================

context("example models")

test_that("definitons for model 1 coincide", {
  data("genePdmp1")
  data("genePoly1")
  
  t <- c(from = 0, to = 10, by = 0.2)
  times(genePdmp1) <- t
  times(genePoly1) <- t
  
  expect_identical(sim(genePdmp1, outSlot = FALSE, seed = 30),
                   sim(genePoly1, outSlot = FALSE, seed = 30))
  })

test_that("definitons for model 2 coincide", {
  data("genePdmp2")
  data("genePoly2")
  
  t <- c(from = 0, to = 10, by = 0.2)
  times(genePdmp2) <- t
  times(genePoly2) <- t
  
  expect_identical(sim(genePdmp2, outSlot = FALSE, seed = 50),
                   sim(genePoly2, outSlot = FALSE, seed = 50))
})

test_that("definitons for model 4 coincide", {
  data("genePdmp4")
  data("genePoly4")
  
  t <- c(from = 0, to = 10, by = 0.2)
  times(genePdmp4) <- t
  times(genePoly4) <- t
  
  expect_identical(sim(genePdmp4, outSlot = FALSE, seed = 40),
                   sim(genePoly4, outSlot = FALSE, seed = 40))
})

test_that("definitons for model 7 coincide", {

  data("genePdmp7")
  data("genePoly7")
  data("toggleSwitch")
  
  t <- c(from = 0, to = 10, by = 0.2)
  times(genePdmp7) <- t
  times(genePoly7) <- t
  times(toggleSwitch) <- t
  
  expect_identical(sim(genePdmp7, outSlot = FALSE, seed = 10),
                   sim(genePoly7, outSlot = FALSE, seed = 10))
  expect_equal(sim(genePdmp7, outSlot = FALSE, seed = 20)[, c("ξA", "ξB")],
            sim(toggleSwitch, outSlot = FALSE, seed = 20)[, c("fA", "fB")],
            check.attributes = FALSE)
})

test_that("definitons for model 8 coincide", {
  data("genePdmp8")
  data("genePoly8")
  
  t <- c(from = 0, to = 10, by = 0.2)
  times(genePdmp8) <- t
  times(genePoly8) <- t
  
  expect_equal(sim(genePdmp8, outSlot = FALSE, seed = 80),
               sim(genePoly8, outSlot = FALSE, seed = 80))
})

test_that("definitons for model Benaim coincide", {
  data("Benaim")
  data("polyBenaim")
  
  t <- c(from = 0, to = 10, by = 0.2)
  times(Benaim) <- t
  times(polyBenaim) <- t
  
  expect_equal(sim(Benaim, outSlot = FALSE, seed = 10),
               sim(polyBenaim, outSlot = FALSE, seed = 10))
})

test_that("definitons for simplePdmp coincide", {
  data("simplePoly")
  data("simplePdmp")
  
  t <- c(from = 0, to = 10, by = 0.2)
  times(simplePdmp) <- t
  times(simplePoly) <- t
  
  expect_identical(sim(simplePoly, outSlot = FALSE, seed = 5),
                   sim(simplePdmp, outSlot = FALSE, seed = 5))
})