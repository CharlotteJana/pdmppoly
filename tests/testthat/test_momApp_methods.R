context("momApp_methods")

data(genePolyF)
data(genePdmpF)
data(genePolyK2)
times(genePolyF) <- c(from = 0, to = 10, by = 1)
times(genePdmpF) <- c(from = 0, to = 10, by = 1)

ma <- momApp(genePolyF, maxorder = 4, 
             closure = c("zero", "normal"), 
             centralize = c(FALSE, TRUE))
ma2 <- momApp(genePolyK2, maxorder = 6, closure = "zero", centralize = TRUE)


test_that("'addSimulation' adds simulations", {

  ma <- addSimulation(ma, multSim = multSim(genePdmpF, seeds = 1:2))
  
  expect_identical(levels(ma$moments$method), 
                   c("normal (central)", "simulation", "zero (raw)"))
  expect_equal(nrow(subset(ma$moments, order == 4)), 33)
})

test_that("'addSimulation' throws errors", {
  expect_error(addSimulation(ma, multSim = multSim(genePolyK2, seeds = 3)))
  expect_error(addSimulation(ma, multSim = multSim(genePolyF, seeds = NULL)))
})

test_that("'print' prints an overview", {
  expect_output(print(ma), 
                regexp = "Moment approximation for moments of order > 4")
  expect_output(print(ma), 
                regexp = "Description: Model F: positive feedback")
  expect_output(print(ma2), 
                regexp = c("no closure     6  100 0.5 7.311301e"))
})

test_that("'print' prints an overview", {
  expect_output(print(ma), 
                regexp = "Moment approximation for moments of order > 4")
  expect_output(print(ma), 
                regexp = "Description: Model F: positive feedback")
  expect_output(print(ma2), 
                regexp = c("no closure     6  100 0.5 7.311301e\\+08"))
})

test_that("'summary' prints a summary", {
  expect_output(summary(ma), 
                regexp = "\\$maxorder \\t 4")
  expect_output(summary(ma), 
                regexp = "\\$moments, order =  2 , method =  zero \\(raw\\)")
  expect_output(summary(ma), 
                regexp = "Mean   :0.4977   Mean   :15409")
})

test_that("'tail' prints tail of every data.frame", {
  expect_output(tail(ma), 
                regexp = "normal \\(central\\)     4   10 0.4982533 398500.3")
  expect_output(tail(ma), 
                regexp = "\\$out\\$zero \\(raw\\)")
  expect_output(tail(ma), 
                regexp = "9 0.5019206 0.4980794 0.5149376 14.17783   9.703366")
})

test_that("'head' prints head of every data.frame", {
  expect_output(head(ma2), 
                regexp = "no closure     1  0.3 0.5 0.6457029 0.08569497")
  expect_output(head(ma2), 
                regexp = "E\\(f1\\*f2\\^5\\|d=0\\)")
  expect_output(head(ma2), 
                regexp = "\\$out\\$no closure")
})
