library(spray)
#------ code to generate the pdmpModel version -----

genePdmpKF <- new("pdmpModel",
                 descr = "Model KF: positive feedback with constant rate",
                 parms = list(b = 0.2, a = 7, k10 = 0.02, k01 = 0.02,
                               m01 = 0.01, m10 = 0.01), 
                 init = c(f = 1, d = 1),
                 discStates = list(d = 0:1),
                 dynfunc = function(t, x, parms) {
                   df <- with(as.list(c(x, parms)), {a*d - b*f})
                   return(c(df, 0))
                 }, 
                 ratefunc = function(t, x, parms) {
                   return(with(as.list(c(x, parms)), 
                               switch(d + 1, k01*f + m01, k10 + m10)))
                 }, 
                 jumpfunc = function(t, x, parms, jtype) {
                   c(x[1], 1 - x[2])
                 }, 
                 times = c(from = 0, to = 100, by = 0.1), 
                 solver = "lsodar")

#------ code to generate the polyPdmpModel version -----

genePolyKF <- new("polyPdmpModel",
                 descr = "Model KF: positive feedback with constant rate (polynomial version)",
                 parms = list(b = 0.2, a = 7, k10 = 0.02, k01 = 0.02,
                               m01 = 0.01, m10 = 0.01), 
                 init = c(f = 1, d = 1), 
                 discStates = list(d = 0:1),
                 dynpolys = quote(list(
                   list(overall = linear(c(-b,a)))
                 )),
                 ratepolys = quote(list(  
                   list(k01*lone(1,2) + m01, k10 + m10)
                 )),
                 jumpfunc = function(t, x, parms, jtype) {
                   c(x[1], 1 - x[2])
                 }, 
                 times = c(from = 0, to = 100, by = 0.1), 
                 solver = "lsodar")

#------- comparison of the models --------------

identical(sim(genePdmpKF, outSlot = FALSE, seed = 20),
          sim(genePolyKF, outSlot = FALSE, seed = 20))
