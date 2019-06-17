#======== todo =================================================================
# v1: verschiedene Ergebnisse von 2b/4b-unimodal für momApp(genePolyF, maxorder = 4, closure = "zero", centralize = TRUE)

context("modalityTest")

data(genePolyF)
init(genePolyF) <- c(f = 0, d = 1)
times(genePolyF) <- c(from = 0, to = 500, by = 1)

ma <- momApp(genePolyF, maxorder = 4, closure = "zero", centralize = TRUE)
sup <- getSupport(genePolyF)

test_that("'modalityTest' allows different inputs for 'lower' and 'upper'", {
  
  res1 <- suppressWarnings(
    modalityTest(ma, lower = sup$lower, upper = sup$upper))
  res2 <- suppressWarnings(
    modalityTest(ma, lower = c(0, 0), upper = sup$upper))
  res3 <- suppressWarnings(
    modalityTest(ma, lower = sup$lower, upper = c(35, 1)))
  res4 <- suppressWarnings(
    modalityTest(ma, lower = c(0, 0), upper = c(35, 1)))
  
  expect_identical(res1, res2)
  expect_identical(res3, res4)
  expect_identical(tail(res1), tail(res3))
  
  moments <- t(subset(ma$moments, time == 10)[, 5:4])
  modality <- is.unimodal(moments, lower = c(0, 0), upper = c(35, 1))
  suppressWarnings(
  expect_equal(dplyr::mutate_if(subset(res4, time == 10), 
                                is.factor,
                                as.character), 
               dplyr::mutate_if(data.frame(time = 10, 
                                           method = "zero (central)",
                                           variable = c("f", "d"), 
                                           modality = modality),
                                is.factor,
                                as.character))
  )
})

test_that("argument 'vars' of 'modalityTest' works as expected", {
  
  res1 <- suppressWarnings(
    modalityTest(ma, lower = sup$lower, upper = sup$upper, vars = "f"))
  res2 <- suppressWarnings(
    modalityTest(ma, lower = 0, upper = 35, vars = "f"))
  
  expect_identical(tail(res1), tail(res2))
  expect_equal(dplyr::mutate_if(subset(res1, time == 10), 
                                is.factor,
                                as.character), 
               dplyr::mutate_if(data.frame(time = 10, 
                                           method = "zero (central)",
                                           variable = "f", 
                                           modality = "not existant"),
                                is.factor,
                                as.character))
})

test_that("'modalityTest' throws errors and warnings", {
  expect_error(modalityTest(ma))
  expect_error(modalityTest(ma, lower = c(0, 0), upper = 35))
  expect_error(modalityTest(ma, lower = c(0, 0), upper = c(35, 1), vars = "f"))
  expect_error(modalityTest(ma, lower = sup$lower[, 1:2], upper = 35, vars = "d"))
})

