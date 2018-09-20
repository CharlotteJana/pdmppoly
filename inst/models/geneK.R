library(spray)
#------ code to generate the pdmpModel version -----

genePdmpK <- new("pdmpModel",
   descr = "Model K: constant activation",
   parms = list(β = 0.005, α = 1, κ10 = 0.01, κ01 = 0.01),
   init = c(ξ = 0, θ = 1),
   discStates = list(θ = 0:1),
   dynfunc = function(t, x, parms) {
     dξ <- with(as.list(c(x, parms)), α*θ - β*ξ)
     return(c(dξ, 0))
   },
   ratefunc = function(t, x, parms) {
     return(with(as.list(c(x, parms)), switch(θ + 1, κ01, κ10)))
   },
   jumpfunc = function(t, x, parms, jtype) {
     c(x[1], 1 - x[2])
   },
   times = c(from = 0, to = 100, by = 0.1),
   solver = "lsodar")

#------ code to generate the polyPdmpModel version -----

genePolyK <- new("polyPdmpModel",
     descr = "Model K: constant activation (polynomial version)",
     parms = list(β = 0.005, α = 1, κ10 = 0.01, κ01 = 0.01),
     init = c(ξ = 0, θ = 1),
     discStates = list(θ = 0:1),
     dynpolys = quote(list(
       list(overall = linear(c(-β, α)))
     )),
     ratepolys = quote(list(
       list(κ01, κ10)
     )),
     jumpfunc = function(t, x, parms, jtype) {
       c(x[1], 1 - x[2])
     },
     times = c(from = 0, to = 100, by = 0.1),
     solver = "lsodar")

#------- comparison of the models --------------

identical(sim(genePdmpK, outSlot = FALSE, seed = 20),
          sim(genePolyK, outSlot = FALSE, seed = 20))
