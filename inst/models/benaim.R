library(spray)
#------ code to generate the pdmpModel version -----

Benaim <- new("pdmpModel",
    descr = "Modell B",
    parms = list(b = 2, β = 1.4),
    init = c(f1 = 10, f2 = 10, d = 1), 
    discStates = list(d = 0:1),
    dynfunc = function(t, x, parms) {
      df <- with(as.list(c(x, parms)), 
                 c(-f1, -f2) + switch(d+1, c(2*b*f2, 0), 
                                       c(0, 2*b*f1)))
      return(c(df, 0))
    }, 
    ratefunc = function(t, x, parms) {
      return(with(as.list(c(x, parms)), β/2))
    }, 
    jumpfunc = function(t, x, parms, jtype) {
      c(x[1:2], 1 - x[3])
    }, 
    times = c(from = 0, to = 100, by = 0.1),
    solver = "lsodar")

#------ code to generate the polyPdmpModel version -----

polyBenaim <- new("polyPdmpModel",
    descr = "Modell B",
    parms = list(b = 2, β = 1.4),
    init = c(f1 = 10, f2 = 10, d = 1), 
    discStates = list(d = 0:1),
    dynpolys = quote(list(
      list(overall = -lone(1,3), specific = list(2*b*lone(2,3), 0)),
      list(overall = -lone(2,3), specific = list(0, 2*b*lone(1,3)))
    )), 
    ratepolys = quote(list(
      list(β/2, β/2)
    )), 
    jumpfunc = function(t, x, parms, jtype) {
      c(x[1:2], 1 - x[3])
    }, 
    times = c(from = 0, to = 100, by = 0.1),
    solver = "lsodar")

#------- comparison of the models --------------

identical(sim(Benaim, outSlot = FALSE, seed = 10),
          sim(polyBenaim, outSlot = FALSE, seed = 10))
