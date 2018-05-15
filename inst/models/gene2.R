#------ code to generate the pdmpModel version -----

genePdmp2 <- new("pdmpModel",
   descr = "Model 2: constant activation with translation",
   parms = list(β1 = 0.025, α1 = 1, κ10 = 0.01, κ01 = 0.01, α2 = 0.5, β2 = 0.02),
   init = c(ξ1 = 0.5, ξ2 = 0, θ = 1), 
   discStates = list(θ = 0:1),
   dynfunc = function(t, x, parms) {
     dξ <- with(as.list(c(x, parms)), c(α1*θ - β1*ξ1, α2*ξ1-β2*ξ2))
     return(c(dξ, 0))
   }, 
   ratefunc = function(t, x, parms) {
     return(with(as.list(c(x, parms)), switch(θ + 1, κ01, κ10)))
   }, 
   jumpfunc = function(t, x, parms, jtype) {
     c(x[1:2], 1 - x[3])
   }, 
   times = c(from = 0, to = 100, by = 0.1), 
   solver = "lsodar")

#------ code to generate the polyPdmpModel version -----

genePoly2 <- new("polyPdmpModel",
   descr = "Model 2: constant activation with translation (polynomial version)",
   parms = list(β1 = 0.025, α1 = 1, κ10 = 0.01, κ01 = 0.01, α2 = 0.5, β2 = 0.02),
   init = c(ξ1 = 0.5, ξ2 = 0, θ = 1), 
   discStates = list(θ = 0:1),
   dynpolys = quote(list(
     list(overall = linear(c(-β1, 0, α1))), #dξ1/dt = -β1ξ1 + α1θ
     list(overall = linear(c(α2, -β2, 0)))  #dξ2/dt = α2ξ1 - β2ξ2
   )),
   ratepolys = quote(list(
     list(κ01,κ10)
   )),
   jumpfunc = function(t, x, parms, jtype) {
     c(x[1:2], 1 - x[3])
   }, 
   times = c(from = 0, to = 100, by = 0.1),
   solver = "lsodar")

#------- comparison of the models --------------

identical(sim(genePdmp2, outSlot = FALSE, seed = 20),
          sim(genePoly2, outSlot = FALSE, seed = 20))
