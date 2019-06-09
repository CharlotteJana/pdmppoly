#======== todo =================================================================
#v2 sollte getSupport für Modell D wirklich einen Fehler zurückgeben?

context("supports")

test_that("getSupport works for models K, F, KF and BF", {

  ### set values
  
  data(genePolyK)
  data(genePdmpF) 
  data(genePolyKF)
  data(genePdmpBF)
  
  a <- 1
  b <- 2
  parms(genePolyK)["a"] <- a
  parms(genePolyK)["b"] <- b
  parms(genePdmpF)["a"] <- a
  parms(genePdmpF)["b"] <- b
  parms(genePolyKF)["a"] <- a
  parms(genePolyKF)["b"] <- b
  parms(genePdmpBF)["a1"] <- a
  parms(genePdmpBF)["a0"] <- 0
  parms(genePdmpBF)["b"] <- b
  
  times <- c(from = 90, to = 100, by = 0.1)
  t <- fromtoby(times)
  times(genePolyK) <- times
  times(genePdmpF) <- times
  times(genePolyKF) <- times
  times(genePdmpBF) <- times
  
  z <- 10
  init(genePolyK)["f"] <- z
  init(genePdmpF)["f"] <- z
  init(genePolyKF)["f"] <- z
  init(genePdmpBF)["f"] <- z
  
  ### create expected result
  
  lower <- data.frame(time = t,
                      f = z*exp(-b*t),
                      d = 0)
  upper <- data.frame(time = t,
                      f = z*exp(-b*t)+a/b*(1-exp(-b*t)),
                      d = 1)
  result <- list(lower = lower, upper = upper)
  
  ### compare results
  
  expect_identical(result, getSupport(genePolyK))
  expect_identical(result, getSupport(genePdmpF))
  expect_identical(result, getSupport(genePolyKF))
  expect_identical(result, getSupport(genePdmpBF))
})

test_that("getSupport throws error for model DF", {
  data(genePolyDF)
  data(Benaim)
  expect_error(getSupport(genePolyDF))
  expect_error(getSupport(Benaim))
})

test_that("getSupport works for model T",{

  data(genePdmpT)
  
  aA <- parms(genePdmpT)[["aA"]]
  bA <- parms(genePdmpT)[["bA"]]
  aB <- parms(genePdmpT)[["aB"]]
  bB <- parms(genePdmpT)[["bB"]]
  
  times <- c(from = 9990, to = 10000, by = 0.1) # we only check the limes
  t <- fromtoby(times)
  times(genePdmpT) <- times
  
  result <- list(lower = data.frame(time = t, fA = 0, fB = 0, d = 1),
                 upper = data.frame(time = t, fA = aA/bA, fB = aB/bB, d = 4))
  
  expect_equal(result, getSupport(genePdmpT))
})

test_that("getSupport works for models K2", {
  
  ### set values
  
  data(genePolyK2)

  a1 <- 2
  b1 <- 3
  a2 <- 4
  b2 <- 5
  parms(genePolyK2)["a1"] <- a1
  parms(genePolyK2)["b1"] <- b1
  parms(genePolyK2)["a2"] <- a2
  parms(genePolyK2)["b2"] <- b2
  
  times <- c(from = 90, to = 100, by = 0.1)
  t <- fromtoby(times)
  times(genePolyK2) <- times
  
  z1 <- 4
  z2 <- 7
  init(genePolyK2)["f1"] <- z1
  init(genePolyK2)["f2"] <- z2
  
  ### create expected result
  
  A <- z2*exp(-b2*t) + (exp(-b1*t)-exp(-b2*t))*(z1*a2)/(b2-b1)
  lower <- data.frame(time = t,
                      f1 = z1*exp(-b1*t),
                      f2 = A,
                      d = 0)
  upper <- data.frame(time = t,
                      f1 = z1*exp(-b1*t)+a1/b1*(1-exp(-b1*t)),
                      f2 = A + (a1*a2)/(b1*b2*(b2-b1))*(b2*(1-exp(-b1*t))-b1*(1-exp(-b2*t))),
                      d = 1)
  result <- list(lower = lower, upper = upper)
  
  ### compare results
  
  expect_equal(result, getSupport(genePolyK2))
})
