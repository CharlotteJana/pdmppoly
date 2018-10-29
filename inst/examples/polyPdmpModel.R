#### a simple PDMP ####

#------ pdmpModel version -----
simplePdmp <- new("pdmpModel",
                  descr = "A simple PDMP",
                  init = c(f = 0, d = 0),
                  discStates = list(d = c(-1, 0, 1)),
                  times = c(from = 0, to = 10, by = 0.01),
                  dynfunc = function(t, x, parms) c(x["d"], 0),
                  ratefunc = function(t, x, parms) c(1+x["d"], 1-x["d"]),
                  jumpfunc = function(t, x, parms, jtype){
                    c(0, switch(jtype, x["d"]-1, x["d"]+1))
                  }
)

#------ polyPdmpModel version -----
simplePoly <- new("polyPdmpModel",
                  descr = "polyModel with two jumptypes",
                  init = c(f = 0, d = 0),
                  times = c(from = 0, to = 10, by = 0.01),
                  discStates = list(d = -1:1),
                  dynpolys = quote(list(
                    list(overall = lone(2,2)) # variant 2
                  )),
                  ratepolys = quote(list(
                    list(0, 1, 2), # first jumptype
                    list(2, 1, 0)  # second jumptype
                  )),
                  jumpfunc = function(t, x, parms, jtype){
                    c(0, switch(jtype, x["d"]-1, x["d"]+1))
                  }
)

#------- comparison of the models --------------
identical(sim(simplePoly, outSlot = FALSE, seed = 5),
          sim(simplePdmp, outSlot = FALSE, seed = 5))

#### the toggleSwitch model ####

#------ pdmpModel version -----
genePdmpT <- new("pdmpModel", 
                 descr = "toggleswitch with two promotors",
                 parms = list(βA = 0.5, βB = 0.5, αA = 2, αB = 4, 
                               κ01A = 0.5, κ10A = 2, κ01B = 0.3, κ10B = 3),
                 init = c(ξA = 0.5, ξB = 0.5, θ = 4),
                 discStates = list(θ = 1:4),
                 dynfunc = function(t, x, parms) {
                   dξ <- with(as.list(c(x, parms)), 
                              c(-βA * ξA, -βB * ξB) + switch(θ, 
                                                              c(0, 0), 
                                                              c(αA, 0), 
                                                              c(0, αB), 
                                                              c(αA, αB)))
                   return(c(dξ, 0))
                 }, 
                 ratefunc = function(t, x, parms) {
                   return(with(as.list(c(x, parms)),
                               c(switch(θ, κ01B, κ01B, κ10B*ξA, κ10B*ξA),
                                 switch(θ, κ01A, κ10A*ξB, κ01A, κ10A*ξB))))
                 }, 
                 jumpfunc = function(t, x, parms, jtype) {
                   c(x[1:2], switch(jtype, 
                                    switch(x[3], 3, 4, 1, 2), 
                                    switch(x[3], 2, 1, 4, 3)))
                 }, 
                 times = c(from = 0, to = 100, by = 0.01), 
                 solver = "lsodar")

#------ polyPdmpModel version -----
library("spray")
genePolyT <- new("polyPdmpModel",
                 descr = "toggleswitch with two promotors (polynomial version)",
                 parms = list(βA = 0.5, βB = 0.5, αA = 2, αB = 4, 
                               κ01A = 0.5, κ10A = 2, κ01B = 0.3, κ10B = 3),
                 init = c(ξA = 0.5, ξB = 0.5, θ = 4), 
                 discStates = list(θ = 1:4),
                 dynpolys = quote(list(
                   list(overall = -βA*lone(1,3), specific = list(0, αA, 0, αA)),
                   list(overall = -βB*lone(2,3), specific = list(0, 0, αB, αB))
                 )), 
                 ratepolys = quote(list(  
                   list(κ01B, κ01B, κ10B*lone(1,3), κ10B*lone(1,3)),
                   list(κ01A, κ10A*lone(2,3), κ01A, κ10A*lone(2,3))
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
all.equal(sim(genePdmpT, outSlot = FALSE, seed = 20)[, c("ξA", "ξB")],
          sim(toggleSwitch, outSlot = FALSE, seed = 20)[, c("fA", "fB")],
          check.attributes = FALSE)
