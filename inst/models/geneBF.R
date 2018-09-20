library(spray)
#------ code to generate the pdmpModel version -----

genePdmpBF <- new("pdmpModel",
                  descr = "Model BF: positive feedback with basal transcription",
                  parms = list(β = 0.2, α0 = 1, α1 = 7, κ10 = 0.02, κ01 = 0.02), 
                  init = c(ξ = 1, θ = 1),
                  discStates = list(θ = 0:1),
                  dynfunc = function(t, x, parms) {
                    dξ <- with(as.list(c(x, parms)), {
                      switch(θ+1, α0 - β*ξ, α1 - β*ξ)
                    })
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

genePolyBF <- new("polyPdmpModel",
                  descr = "Model BF: positive feedback with basal transcription (polynomial version)",
                  parms = list(β = 0.2, α0 =1, α1 = 7, κ10 = 0.02, κ01 = 0.02), 
                  init = c(ξ = 1, θ = 1), 
                  discStates = list(θ = 0:1),
                  dynpolys = quote(list(
                    list(overall = -β*lone(1,2),
                         specific = list(α0, α1))
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

identical(sim(genePdmpBF, outSlot = FALSE, seed = 20),
          sim(genePolyBF, outSlot = FALSE, seed = 20))
