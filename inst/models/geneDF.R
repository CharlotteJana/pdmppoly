
#------ code to generate the pdmpModel version -----

genePdmpDF <- new("pdmpModel", 
     descr = "Model DF: dimers + positive feedback",
     parms = list(b = 1, a = 1, k01 = 1, k10 = 1, m21 = 1, m12 = 1),
     init = c(f = 1, fd = 0.5, d = 1), 
     discStates = list(d = 0:1),
     dynfunc = function(t, x, parms) {
       df <- with(as.list(c(x, parms)), 
                  c(-2*m21*f^2 + 2*m12*fd - b*f + d*a, 
                    m21*f^2 - m12*fd))
       return(c(df, 0))
     }, 
     ratefunc = function(t, x, parms) {
       return(with(as.list(c(x, parms)),
                   c(switch(d+1, k01*fd, k10))))
     }, 
     jumpfunc = function(t, x, parms, jtype) {
       c(x[1:2], 1 - x[3])
     }, 
     times = c(from = 0, to = 100, by = 0.1), 
     solver = "lsodar")

#------ code to generate the polyPdmpModel version -----

library("spray")
genePolyDF <- new("polyPdmpModel", 
    descr = "Model DF: dimers + positive feedback (polynomial version)",
    parms = list(b = 1, a = 1, k01 = 1, k10 = 1, m21 = 1, m12 = 1), 
    init = c(f = 1, fd = 0.5, d = 1), 
    discStates = list(d = 0:1),
    dynpolys = quote(list(
      list(overall = linear(c(-2*m21, 0, 0), 2) + linear(c(-b, 2*m12, a))),
      list(overall = linear(c(m21, 0, 0), 2) - m12*lone(2,3))
    )),
    ratepolys = quote(list(
      list(k01*lone(2,3), k10)
    )),
    jumpfunc = function(t, x, parms, jtype){
      c(x[1:2], 1 - x[3])
    }, 
    times = c(from = 0, to = 100, by = 0.1), 
    solver = "lsodar")

#------- comparison of the models --------------

all.equal(sim(genePdmpDF, outSlot = FALSE, seed = 12),
          sim(genePolyDF, outSlot = FALSE, seed = 12))

