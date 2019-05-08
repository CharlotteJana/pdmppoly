#======== todo =================================================================

context("supports")

test_that("getSupport returns a data.frame", {
  data(genePolyK)
  a <- 1
  b <- 0.005
  times(genePolyK) <- c(from = 9900, to = 10000, by = 0.1)
  parms(genePolyK) <- c(b = b, a = a, k10 = 0.01, k01 = 0.01)
  s1 <- getSupport(genePolyK)
  s2 <- data.frame(times = 10000.0, lower = 0, upper = a/b)
  expect_equal(s1[nrow(s1), ], s2, check.attributes = FALSE)
  
})

test_that("results of getSupport are identical for pdmpModel and polyPdmpModel",{
  data(genePolyT)
  data(genePdmpT)
  
  s1 <- getSupport(genePolyT)
  s2 <- getSupport(genePdmpT)
  
  expect_identical(colnames(s1), c("times","lower1","lower2","upper1","upper2"))
  expect_identical(s1, s2)
})
