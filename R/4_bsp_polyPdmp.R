####################################################
########    Examples for polynomial PDMPs   ########
####################################################

#' @include polypdmp_class.R polypdmp_accessors.R

#### model1: constant activation (only transcription) #####
polyModel1 <- polyPdmpModel(
              descr = "Model 1: constant activation (only transcription)", 
              discStates = list(θ = 0:1),
              parms = list(β = 1, α = 2, κ10 = 0.5, κ01 = 3),
              init = c(ξ = 0.5, θ = 1), 
              dynpolys = quote(list(
                            list(overall = linear(c(-β,α))) # one continuous variable, dξ/dt = -β*ξ + α*θ 
                            # alternative: list(-β*lone(1,2), -β*lone(1,2)+α)            
                )), 
              ratepolys = quote(list(  
                            list(κ01, κ10) # one jumptype, constant rates κ01 and κ10
                          )), 
              jumpfunc = function(t, x, parms, jtype) {
                c(x[1], 1 - x[2])
              }, 
              times = c(from = 0, to = 10, by = 0.01), 
              solver = "lsodar")


#### model7: toggleswitch with two promotors ####
#   θ=1 A blocks B, B blocks A 
#   θ=2 A blocks B, B not A 
#   θ=3 A not    B, B blocks A
#   θ=4 nothing blocked
#   jtype=1: promotor B changes (A binds/unbinds) 
#   jtype=2: promotor A changes (B binds/unbinds)
polyModel7 <- polyPdmpModel(
              descr = "toggleswitch with one discrete Variable",
              discStates = list(θ = 1:4),
              parms = list(βA = 0.5, βB = 1,   αA = 2, αB = 4, 
                            κ01A = 0.5, κ10A = 2, κ01B = 1/3, κ10B = 3), 
              init = c(ξA = 0.5, ξB = 0.5, θ = 4), 
              dynpolys = quote(list(
                            list(overall = -βA*lone(1,3), 
                                 specific = list(0, αA, 0, αA)), # dξA/dt for θ = 1,...,4
                            list(overall = -βB*lone(2,3), 
                                 specific = list(0, 0, αB, αB))  # dξB/dt for θ = 1,...,4
                          )), 
              ratepolys = quote(list(  
                            list(κ01A, κ01A, κ10A*lone(2,3), κ10A*lone(2,3)), # first jumptype
                            list(κ01B, κ10B*lone(1,3), κ01B, κ10B*lone(1,3))  # second jumptype
                          )),
              jumpfunc = function(t, x, parms, jtype) {
                c(x[1:2], switch(jtype, 
                                 switch(x[3], 3, 4, 1, 2), 
                                 switch(x[3], 2, 1, 4, 3)))
              }, 
              times = c(from = 0, to = 10, by = 0.01), 
              solver = "lsodar")
