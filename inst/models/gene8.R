#------ code to generate the pdmpModel version -----

genePdmp8 <- new("pdmpModel", 
     descr = "Model 8: dimers + positive feedback",
     parms = list(β = 1, α = 1, κ01 = 1, κ10 = 1, γ21 = 1, γ12 = 1),
     init = c(ξ = 1, ξd = 0.5, θ = 1), 
     discStates = list(θ = 0:1),
     dynfunc = function(t, x, parms) {
       dξ <- with(as.list(c(x, parms)), 
                  c(-2*γ21*ξ^2 + 2*γ12*ξd - β*ξ + θ*α, 
                    γ21*ξ^2 - γ12*ξd))
       return(c(dξ, 0))
     }, 
     ratefunc = function(t, x, parms) {
       return(with(as.list(c(x, parms)),
                   c(switch(θ+1, κ01*ξd, κ10))))
     }, 
     jumpfunc = function(t, x, parms, jtype) {
       c(x[1:2], 1 - x[3])
     }, 
     times = c(from = 0, to = 100, by = 0.01), 
     solver = "lsodar")

#------ code to generate the polyPdmpModel version -----

library("spray")
genePoly8 <- new("polyPdmpModel", 
    descr = "Model 8: dimers + positive feedback (polynomial version)",
    parms = list(β = 1, α = 1, κ01 = 1, κ10 = 1, γ21 = 1, γ12 = 1), 
    init = c(ξ = 1, ξd = 0.5, θ = 1), 
    discStates = list(θ = 0:1),
    dynpolys = quote(list(
      list(overall = linear(c(-2*γ21, 0, 0), 2) + linear(c(-β, 2*γ12, α))),
      list(overall = linear(c(γ21, 0, 0), 2) - γ12*lone(2,3))
    )),
    ratepolys = quote(list(
      list(κ01*lone(2,3), κ10)
    )),
    jumpfunc = function(t, x, parms, jtype){
      c(x[1:2], 1 - x[3])
    }, 
    times = c(from = 0, to = 100, by = 0.01), 
    solver = "lsodar")

#------- comparison of the models --------------

all.equal(sim(genePdmp8, outSlot = FALSE, seed = 12),
          sim(genePoly8, outSlot = FALSE, seed = 12))
