library(spray)
#------ code to generate the pdmpModel version -----

genePdmpK <- new("pdmpModel",
   descr = "Model K: constant activation",
   parms = list(b = 0.005, a = 1, k10 = 0.01, k01 = 0.01),
   init = c(f = 0, d = 1),
   discStates = list(d = 0:1),
   dynfunc = function(t, x, parms) {
     df <- with(as.list(c(x, parms)), a*d - b*f)
     return(c(df, 0))
   },
   ratefunc = function(t, x, parms) {
     return(with(as.list(c(x, parms)), switch(d + 1, k01, k10)))
   },
   jumpfunc = function(t, x, parms, jtype) {
     c(x[1], 1 - x[2])
   },
   times = c(from = 0, to = 100, by = 0.1),
   solver = "lsodar")

#------ code to generate the polyPdmpModel version -----

genePolyK <- new("polyPdmpModel",
     descr = "Model K: constant activation (polynomial version)",
     parms = list(b = 0.005, a = 1, k10 = 0.01, k01 = 0.01),
     init = c(f = 0, d = 1),
     discStates = list(d = 0:1),
     dynpolys = quote(list(
       list(overall = linear(c(-b, a)))
     )),
     ratepolys = quote(list(
       list(k01, k10)
     )),
     jumpfunc = function(t, x, parms, jtype) {
       c(x[1], 1 - x[2])
     },
     times = c(from = 0, to = 100, by = 0.1),
     solver = "lsodar")

#------- comparison of the models --------------

identical(sim(genePdmpK, outSlot = FALSE, seed = 20),
          sim(genePolyK, outSlot = FALSE, seed = 20))
