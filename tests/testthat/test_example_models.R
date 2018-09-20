#======== todo =================================================================

context("example models")

test_that("definitons for model K coincide", {
  data("genePdmpK")
  data("genePolyK")
  
  t <- c(from = 0, to = 10, by = 0.2)
  times(genePdmpK) <- t
  times(genePolyK) <- t
  
  expect_identical(sim(genePdmpK, outSlot = FALSE, seed = 30),
                   sim(genePolyK, outSlot = FALSE, seed = 30))
  })

test_that("definitons for model K2 coincide", {
  data("genePdmpK2")
  data("genePolyK2")
  
  t <- c(from = 0, to = 10, by = 0.2)
  times(genePdmpK2) <- t
  times(genePolyK2) <- t
  
  expect_identical(sim(genePdmpK2, outSlot = FALSE, seed = 50),
                   sim(genePolyK2, outSlot = FALSE, seed = 50))
})

test_that("definitons for model F coincide", {
  data("genePdmpF")
  data("genePolyF")
  
  t <- c(from = 0, to = 10, by = 0.2)
  times(genePdmpF) <- t
  times(genePolyF) <- t
  
  expect_identical(sim(genePdmpF, outSlot = FALSE, seed = 40),
                   sim(genePolyF, outSlot = FALSE, seed = 40))
})

test_that("definitons for model KF coincide", {
  data("genePdmpKF")
  data("genePolyKF")
  
  t <- c(from = 0, to = 10, by = 0.2)
  times(genePdmpKF) <- t
  times(genePolyKF) <- t
  
  expect_identical(sim(genePdmpKF, outSlot = FALSE, seed = 10),
                   sim(genePolyKF, outSlot = FALSE, seed = 10))
})

test_that("definitons for model T coincide", {

  data("genePdmpT")
  data("genePolyT")
  data("toggleSwitch")
  
  t <- c(from = 0, to = 10, by = 0.2)
  times(genePdmpT) <- t
  times(genePolyT) <- t
  times(toggleSwitch) <- t
  
  expect_identical(sim(genePdmpT, outSlot = FALSE, seed = 10),
                   sim(genePolyT, outSlot = FALSE, seed = 10))
  expect_equal(sim(genePdmpT, outSlot = FALSE, seed = 20)[, c("ξA", "ξB")],
            sim(toggleSwitch, outSlot = FALSE, seed = 20)[, c("fA", "fB")],
            check.attributes = FALSE)
})

test_that("definitons for model DF coincide", {
  data("genePdmpDF")
  data("genePolyDF")
  
  t <- c(from = 0, to = 10, by = 0.2)
  times(genePdmpDF) <- t
  times(genePolyDF) <- t
  
  expect_equal(sim(genePdmpDF, outSlot = FALSE, seed = 80),
               sim(genePolyDF, outSlot = FALSE, seed = 80))
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