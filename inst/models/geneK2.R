library(spray)
#------ code to generate the pdmpModel version -----

genePdmpK2 <- new("pdmpModel",
   descr = "Model K2: constant activation with translation",
   parms = list(b1 = 0.025, a1 = 1, k10 = 0.01, k01 = 0.01, a2 = 0.5, b2 = 0.02),
   init = c(f1 = 0.5, f2 = 0, d = 1), 
   discStates = list(d = 0:1),
   dynfunc = function(t, x, parms) {
     df <- with(as.list(c(x, parms)), c(a1*d - b1*f1, a2*f1-b2*f2))
     return(c(df, 0))
   }, 
   ratefunc = function(t, x, parms) {
     return(with(as.list(c(x, parms)), switch(d + 1, k01, k10)))
   }, 
   jumpfunc = function(t, x, parms, jtype) {
     c(x[1:2], 1 - x[3])
   }, 
   times = c(from = 0, to = 100, by = 0.1), 
   solver = "lsodar")

#------ code to generate the polyPdmpModel version -----

genePolyK2 <- new("polyPdmpModel",
   descr = "Model K2: constant activation with translation (polynomial version)",
   parms = list(b1 = 0.025, a1 = 1, k10 = 0.01, k01 = 0.01, a2 = 0.5, b2 = 0.02),
   init = c(f1 = 0.5, f2 = 0, d = 1), 
   discStates = list(d = 0:1),
   dynpolys = quote(list(
     list(overall = linear(c(-b1, 0, a1))), #df1/dt = -b1f1 + a1d
     list(overall = linear(c(a2, -b2, 0)))  #df2/dt = a2f1 - b2f2
   )),
   ratepolys = quote(list(
     list(k01,k10)
   )),
   jumpfunc = function(t, x, parms, jtype) {
     c(x[1:2], 1 - x[3])
   }, 
   times = c(from = 0, to = 100, by = 0.1),
   solver = "lsodar")

#------- comparison of the models --------------

identical(sim(genePdmpK2, outSlot = FALSE, seed = 20),
          sim(genePolyK2, outSlot = FALSE, seed = 20))
