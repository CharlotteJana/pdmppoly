#======== todo =================================================================

context("polyPdmpModel")

#========== simplePdmp ================

test_that("polyPdmpModel can be defined as pdmpModels - two jumptypes", {
  library(spray)
  simplePdmp <- pdmpModel(
    descr = "Model with two jumptypes",
    init = c(f = 0, d = 0),
    times = c(from = 0, to = 10, by = 0.1),
    discStates = list(d = -1:1),
    dynfunc = function(t, x, parms) c(x["d"], 0),
    ratefunc = function(t, x, parms) c(1+x["d"], 1-x["d"]),
    jumpfunc = function(t, x, parms, jtype){
      c(0, switch(jtype, x["d"]-1, x["d"]+1))
    }
  )
  
  simplePoly <- polyPdmpModel(
    descr = "polyModel with two jumptypes",
    init = c(f = 0, d = 0),
    times = c(from = 0, to = 10, by = 0.1),
    discStates = list(d = -1:1),
    dynpolys = quote(list(
      list(overall = lone(2,2)) # variant 2
      )),
    ratepolys = quote(list(
      list(0, 1, 2), # first jumptype
      list(2, 1, 0) # second jumptype
    )),
    jumpfunc = function(t, x, parms, jtype){
      c(0, switch(jtype, x["d"]-1, x["d"]+1))
    }
  )
 
  expect_identical(sim(simplePdmp, seed = 5, outSlot = FALSE), 
                   sim(simplePoly, seed = 5, outSlot = FALSE))
  
  dynpolys(simplePoly) <- quote(list(list(-1, 0, 1))) # variant 1
  expect_identical(sim(simplePdmp, seed = 8, outSlot = FALSE), 
                   sim(simplePoly, seed = 8, outSlot = FALSE))
})

#========== two continous variables ================

test_that("polyPdmpModel can be defined as pdmpModel - two cont. variables", {
  library(spray)
  modelPdmp2 <- pdmpModel(
    descr = "Model with two continous variables",
    init = c(f1 = 0.5, f2 = 0, d = 1), 
    times = c(from = 0, to = 10, by = 0.1),
    discStates = list(d = 0:1),
    parms = list(a = 1),
    dynfunc = function(t, x, parms){
      with(as.list(c(x, parms)), c(-3*f1 + a*d, a*f1-2*f2, 0))
    },
    ratefunc = function(t, x, parms) 2*(x["d"]+1),
    jumpfunc = function(t, x, parms, jtype){
      c(x[1:2], 1 - x[3])
    }
  )
  
  modelPoly2 <- new("polyPdmpModel", 
    descr = "PolyModel with two continous variables",
    init = c(f1 = 0.5, f2 = 0, d = 1), 
    times = c(from = 0, to = 10, by = 0.1),
    discStates = list(d = 0:1),
    parms = list(a = 1), 
    dynpolys = quote(list( # variant 2
      list(overall = linear(c(-3, 0, a))), #df1/dt = -3f1 + ad
      list(overall = linear(c(a, -2, 0)))  #df2/dt = af1 - 2f2
    )),
    ratepolys = quote(list(
      list(2, 4)
    )),
    jumpfunc = function(t, x, parms, jtype) {
      c(x[1:2], 1 - x[3])
    }
  )
  
  expect_identical(sim(modelPdmp2, seed = 10, outSlot = FALSE), 
                   sim(modelPoly2, seed = 10, outSlot = FALSE))
  
  dynpolys(modelPoly2) <- quote(list( # variant 1
                            list(linear(c(-3, 0, 0)), linear(c(-3, 0, a))),
                            list(linear(c(a, -2, 0)), linear(c(a, -2, 0)))
                            ))
  expect_identical(sim(modelPdmp2, seed = 3, outSlot = FALSE), 
                   sim(modelPoly2, seed = 3, outSlot = FALSE))
})

#==============================

test_that("order of variables in init doesn't matter for sim", {
  data(genePoly7)
  times(genePoly7) <- c(from = 0, to = 10, by = 0.1)
  init(genePoly7) <-c(ξA = 0.5, ξB = 0.5, θ = 4)
  sim1 <- sim(genePoly7, seed = 2, outSlot = FALSE)
  
  model <- genePoly7 <- new("polyPdmpModel",
                            descr = "Model 7 with different order of variables",
                            parms = parms(genePoly7),
                            init = c(ξA = 0.5, θ = 4, ξB = 0.5), 
                            discStates = list(θ = 1:4),
                            dynpolys = quote(list(
                              list(overall = -βA*lone(1,3), specific = list(0, αA, 0, αA)),
                              list(overall = -βB*lone(3,3), specific = list(0, 0, αB, αB))
                            )), 
                            ratepolys = quote(list(  
                              list(κ01B, κ01B, κ10B*lone(1,3), κ10B*lone(1,3)),
                              list(κ01A, κ10A*lone(3,3), κ01A, κ10A*lone(3,3))
                            )),
                            jumpfunc = function(t, x, parms, jtype) {
                              c(x[1], switch(jtype, 
                                               switch(x[2], 3, 4, 1, 2), 
                                               switch(x[2], 2, 1, 4, 3)), x[3])
                            }, 
                            times = times(genePoly7), 
                            solver = "lsodar")
  
  sim2 <- sim(model, seed = 2, outSlot = FALSE)
  expect_identical(sim1[, "ξA"], sim2[, "ξA"])
  expect_identical(sim1[, "θ"], sim2[, "θ"])
})
  