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
