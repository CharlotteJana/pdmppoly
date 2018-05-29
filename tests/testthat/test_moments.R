context("calculations")
setwd(system.file("testdata", package = "pdmpsim", mustWork = TRUE))
ms <- readRDS("test_MultSim.rda")
msCsv <- loadMultSimCsv("test_MultSimCsv.rda")

#-------------- tests ----------------

test_that("method 'mean' and 'moments(..., order = 1)' give equal results", {
  if(identical(find('colmoment'),"package:LaF")){
    means <- mean(ms)
    meansCsv <- mean(msCsv)
    moments <- moments(ms, order = 1)
    momentsCsv <- moments(msCsv, order = 1)
    
    expect_equal(means, moments)
    expect_equal(means, meansCsv, tolerance = 1e-07)
    expect_equal(moments, momentsCsv, tolerance = 1e-07)
    # 'm.' and 'm.Csv' are not exactly equal because 'multSimCsv' only stores 7 digits.
  }
})

test_that("method 'moments' gives equal results for 'multSim' and 'multSimCsv'", {
  if(identical(find('colmoment'),"package:LaF")){
    m2 <- moments(ms, order = 2)
    m2Csv <- moments(msCsv, order = 2)
    expect_equal(m2, m2Csv, tolerance = 1e-07)
    
    m10 <- moments(ms, order = 5)
    m10Csv <- moments(msCsv, order = 5)
    expect_equal(m10, m10Csv, tolerance = 1e-07)
    # 'm.' and 'm.Csv' are not exactly equal because 'multSimCsv' only stores 7 digits.
  }
})