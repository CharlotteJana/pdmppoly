#------ code to generate this model -----

simplePoly <- polyPdmpModel(
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

data("simplePoly")
data("simplePdmp")
identical(sim(simplePoly, outSlot = FALSE, seed = 5),
          sim(simplePdmp, outSlot = FALSE, seed = 5))
