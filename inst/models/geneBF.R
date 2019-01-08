library(spray)
#------ code to generate the pdmpModel version -----

genePdmpBF <- new("pdmpModel",
                  descr = "Model BF: positive feedback with basal transcription",
                  parms = list(b = 0.2, a0 = 1, a1 = 7, k10 = 0.02, k01 = 0.02), 
                  init = c(f = 1, d = 1),
                  discStates = list(d = 0:1),
                  dynfunc = function(t, x, parms) {
                    df <- with(as.list(c(x, parms)), {
                      switch(d+1, a0 - b*f, a1 - b*f)
                    })
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

genePolyBF <- new("polyPdmpModel",
                  descr = "Model BF: positive feedback with basal transcription (polynomial version)",
                  parms = list(b = 0.2, a0 =1, a1 = 7, k10 = 0.02, k01 = 0.02), 
                  init = c(f = 1, d = 1), 
                  discStates = list(d = 0:1),
                  dynpolys = quote(list(
                    list(overall = -b*lone(1,2),
                         specific = list(a0, a1))
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

identical(sim(genePdmpBF, outSlot = FALSE, seed = 20),
          sim(genePolyBF, outSlot = FALSE, seed = 20))
