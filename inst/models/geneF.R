library(spray)
#------ code to generate the pdmpModel version -----

genePdmpF <- new("pdmpModel",
   descr = "Model F: positive feedback",
   parms = list(β = 0.2, α = 7, κ10 = 0.02, κ01 = 0.02), 
   init = c(ξ = 1, θ = 1),
   discStates = list(θ = 0:1),
   dynfunc = function(t, x, parms) {
     dξ <- with(as.list(c(x, parms)), {α*θ - β*ξ})
     return(c(dξ, 0))
   }, 
   ratefunc = function(t, x, parms) {
     return(with(as.list(c(x, parms)), switch(θ + 1, κ01*ξ, κ10)))
   }, 
   jumpfunc = function(t, x, parms, jtype) {
     c(x[1], 1 - x[2])
   }, 
   times = c(from = 0, to = 100, by = 0.1), 
   solver = "lsodar")

#------ code to generate the polyPdmpModel version -----

genePolyF <- new("polyPdmpModel",
     descr = "Model F: positive feedback (polynomial version)",
     parms = list(β = 0.2, α = 7, κ10 = 0.02, κ01 = 0.02), 
     init = c(ξ = 1, θ = 1), 
     discStates = list(θ = 0:1),
     dynpolys = quote(list(
       list(overall = linear(c(-β,α)))
     )),
     ratepolys = quote(list(  
       list(κ01*lone(1,2), κ10)
     )),
     jumpfunc = function(t, x, parms, jtype) {
       c(x[1], 1 - x[2])
     }, 
     times = c(from = 0, to = 100, by = 0.1), 
     solver = "lsodar")

#------- comparison of the models --------------

identical(sim(genePdmpF, outSlot = FALSE, seed = 40),
          sim(genePolyF, outSlot = FALSE, seed = 40))
