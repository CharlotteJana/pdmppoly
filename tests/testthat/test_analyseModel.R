context("analyseModel")

tmp <- tempfile()

data(genePolyK2)
data(genePdmpK2)
descr(genePolyK2) <- "testK2"
descr(genePdmpK2) <- "testK2"

test_that("'analyseModel' creates correct files", {
  
  # sim = TRUE
  analyseModel(polyModel = genePolyK2, model = genePdmpK2,
               seeds = 1:3, dir = tmp, momentorder = NULL,
               sim = TRUE, plot = FALSE, modality = FALSE,
               momApp = FALSE, statistics = FALSE)
  
  expect_identical(list.files(tmp),
                   c("testK2__multSimData.rda", 
                     "testK2__simulations.rda"))
  
  # statistics = TRUE
  
  analyseModel(polyModel = genePolyK2, model = genePdmpK2,
               seeds = 1:3, dir = tmp, momentorder = NULL,
               sim = TRUE, plot = FALSE, modality = FALSE,
               momApp = FALSE, statistics = TRUE,
               funs = c("min", "max"))
  
  expect_identical(list.files(tmp),
                   c("testK2__multSimData.rda", 
                     "testK2__simulations.rda",
                     "testK2__statistics.rda"))
  
  # momApp = TRUE
  analyseModel(polyModel = genePolyK2, model = genePdmpK2,
               seeds = 1:3, dir = tmp, momentorder = c(4, 5),
               sim = FALSE, plot = FALSE, modality = FALSE,
               momApp = TRUE, statistics = FALSE)
  
  expect_identical(list.files(tmp),
                   c("testK2__moments_order<=4.rda",
                     "testK2__moments_order<=5.rda",
                     "testK2__multSimData.rda", 
                     "testK2__simulations.rda",
                     "testK2__statistics.rda"))
  
  # modality = TRUE
  suppressWarnings(
    analyseModel(polyModel = genePolyK2, model = genePdmpK2,
                 seeds = 1:3, dir = tmp, momentorder = c(4, 5),
                 sim = FALSE, plot = FALSE, statistics = FALSE,
                 modality = TRUE, momApp = FALSE,
                 lower = NULL, upper = getSupport(genePdmpK2)$upper)
  )
  
  expect_identical(list.files(tmp),
                   c("testK2__modality_order<=4.rda",
                     "testK2__modality_order<=5.rda",
                     "testK2__moments_order<=4.rda",
                     "testK2__moments_order<=5.rda",
                     "testK2__multSimData.rda", 
                     "testK2__simulations.rda",
                     "testK2__statistics.rda"))
  
  # plot = TRUE
  analyseModel(polyModel = genePolyK2, model = genePdmpK2,
               seeds = 1:3, dir = tmp, momentorder = c(4, 5),
               sim = FALSE, plot = TRUE, modality = FALSE,
               statistics = FALSE, momApp = FALSE, 
               funs = c("min", "max"))
  
  expect_identical(list.files(tmp),
                   c("testK2__boxplot.png",
                     "testK2__densities.png",
                     "testK2__histogram.png",
                     "testK2__modality_order<=4.png",
                     "testK2__modality_order<=4.rda",
                     "testK2__modality_order<=5.png",
                     "testK2__modality_order<=5.rda",
                     "testK2__moments_order<=4.png",
                     "testK2__moments_order<=4.rda",
                     "testK2__moments_order<=5.png",
                     "testK2__moments_order<=5.rda",
                     "testK2__multSimData.rda", 
                     "testK2__plot.png",
                     "testK2__simulations.rda",
                     "testK2__singleSimulations.png",
                     "testK2__statistics_f1.png",
                     "testK2__statistics_f2.png",
                     "testK2__statistics.rda",
                     "testK2__violins.png"))
})

test_that("'analyseModel' throws errors and warnings", {
  
  # change in momentorder without new calculation
  expect_warning(analyseModel(polyModel = genePolyK2, model = genePdmpK2,
                              seeds = 1:3, dir = tmp, momentorder = 10,
                              sim = FALSE, plot = FALSE, modality = FALSE,
                              lower = NULL, upper = getSupport(genePdmpK2)$upper,
                              momApp = FALSE, statistics = FALSE),
                 regexp = "*moments_order<=10*")
  
  # change in seeds without new simulation
  expect_error(analyseModel(polyModel = genePolyK2, model = genePdmpK2,
                            seeds = 1:5, dir = tmp, momentorder = 4,
                            sim = FALSE, plot = FALSE, modality = FALSE,
                            momApp = FALSE, statistics = FALSE))
  
  # change of model parameters
  parms(genePolyK2)[["b2"]] <- parms(genePolyK2)[["b2"]] + 1
  expect_error(analyseModel(polyModel = genePolyK2, model = genePdmpK2, 
                            seeds = 1:3, dir = tmp, momentorder = 4))
  parms(genePdmpK2)[["b2"]] <- parms(genePolyK2)[["b2"]]
  expect_error(analyseModel(polyModel = genePolyK2, model = genePdmpK2, 
                            seeds = 1:3, dir = tmp, momentorder = 4,
                            sim = FALSE))
})

unlink(tmp, recursive = TRUE) # remove all files
