library(spray)
#------ code to generate the pdmpModel version -----

genePdmpF <- new("pdmpModel",
   descr = "Model F: positive feedback",
   parms = list(b = 0.2, a = 7, k10 = 0.02, k01 = 0.02), 
   init = c(f = 1, d = 1),
   discStates = list(d = 0:1),
   dynfunc = function(t, x, parms) {
     df <- with(as.list(c(x, parms)), {a*d - b*f})
     return(c(df, 0))
   }, 
   ratefunc = function(t, x, parms) {
     return(with(as.list(c(x, parms)), switch(d + 1, k01*f, k10)))
   }, 
   jumpfunc = function(t, x, parms, jtype) {
     c(x[1], 1 - x[2])
   }, 
   times = c(from = 0, to = 100, by = 0.1), 
   solver = "lsodar")

#------ code to generate the polyPdmpModel version -----

genePolyF <- new("polyPdmpModel",
     descr = "Model F: positive feedback (polynomial version)",
     parms = list(b = 0.2, a = 7, k10 = 0.02, k01 = 0.02), 
     init = c(f = 1, d = 1), 
     discStates = list(d = 0:1),
     dynpolys = quote(list(
       list(overall = linear(c(-b,a)))
     )),
     ratepolys = quote(list(  
       list(k01*lone(1,2), k10)
     )),
     jumpfunc = function(t, x, parms, jtype) {
       c(x[1], 1 - x[2])
     }, 
     times = c(from = 0, to = 100, by = 0.1), 
     solver = "lsodar")

#------- comparison of the models --------------

identical(sim(genePdmpF, outSlot = FALSE, seed = 40),
          sim(genePolyF, outSlot = FALSE, seed = 40))
