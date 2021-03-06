library(spray)
#------ code to generate the pdmpModel version -----

genePdmpT <- new("pdmpModel", 
   descr = "Model T: toggleswitch with two promotors",
   parms = list(bA = 0.5, bB = 0.5, aA = 2, aB = 4, 
                 k01A = 0.5, k10A = 2, k01B = 0.3, k10B = 3),
   init = c(fA = 0.5, fB = 0.5, d = 4),
   discStates = list(d = 1:4),
   dynfunc = function(t, x, parms) {
     df <- with(as.list(c(x, parms)), 
                c(-bA * fA, -bB * fB) + switch(d, 
                                               c(0, 0), 
                                               c(aA, 0), 
                                               c(0, aB), 
                                               c(aA, aB)))
     return(c(df, 0))
   }, 
   ratefunc = function(t, x, parms) {
     return(with(as.list(c(x, parms)),
                 c(switch(d, k01B, k01B, k10B*fA, k10B*fA),
                   switch(d, k01A, k10A*fB, k01A, k10A*fB))))
   }, 
   jumpfunc = function(t, x, parms, jtype) {
     c(x[1:2], switch(jtype, 
                      switch(x[3], 3, 4, 1, 2), 
                      switch(x[3], 2, 1, 4, 3)))
   }, 
   times = c(from = 0, to = 100, by = 0.01), 
   solver = "lsodar")

#------ code to generate the polyPdmpModel version -----

library("spray")
genePolyT <- new("polyPdmpModel",
  descr = "Model T: toggleswitch with two promotors (polynomial version)",
  parms = list(bA = 0.5, bB = 0.5, aA = 2, aB = 4, 
               k01A = 0.5, k10A = 2, k01B = 0.3, k10B = 3),
  init = c(fA = 0.5, fB = 0.5, d = 4), 
  discStates = list(d = 1:4),
  dynpolys = quote(list(
    list(overall = -bA*lone(1,3), specific = list(0, aA, 0, aA)),
    list(overall = -bB*lone(2,3), specific = list(0, 0, aB, aB))
  )), 
  ratepolys = quote(list(  
    list(k01B, k01B, k10B*lone(1,3), k10B*lone(1,3)),
    list(k01A, k10A*lone(2,3), k01A, k10A*lone(2,3))
  )),
  jumpfunc = function(t, x, parms, jtype) {
    c(x[1:2], switch(jtype, 
                     switch(x[3], 3, 4, 1, 2), 
                     switch(x[3], 2, 1, 4, 3)))
  }, 
  times = c(from = 0, to = 100, by = 0.01), 
  solver = "lsodar")

#------- comparison of the models --------------

identical(sim(genePdmpT, outSlot = FALSE, seed = 10),
          sim(genePolyT, outSlot = FALSE, seed = 10))

data("toggleSwitch")
all.equal(sim(genePdmpT, outSlot = FALSE, seed = 20)[, c("fA", "fB")],
          sim(toggleSwitch, outSlot = FALSE, seed = 20)[, c("fA", "fB")],
          check.attributes = FALSE)
