##### polyPdmpModel ########

setClass("polyPdmpModel", 
         slots = list(discDomain = "integer", dynpolys = "call", ratepolys = "call", dynsprays = "list", ratesprays = "list"), 
         contains = "pdmpModel")

polyPdmpModel <- function(obj = NULL, .cache = new.env(),
                          discDomain, dynpolys, ratepolys, jumpfunc, parms, init, 
                          dynfunc = NULL, ratefunc = NULL, out = NULL, label = character(0), descr = character(0),
                          dynsprays = NULL, ratesprays = NULL,
                          times = c(from = 0, to = 10, by = 1), 
                          solver = "lsodar", ...){
  # see function "initialize" for detailed code
  obj <- new(Class = "polyPdmpModel", discDomain = discDomain, dynpolys = dynpolys, dynfunc = dynfunc,
             jumpfunc = jumpfunc, ratepolys = ratepolys, ratefunc = ratefunc,
             dynsprays = dynsprays, ratesprays = ratesprays, label = label, descr = descr,
             parms = parms, out = out, init = init, times = times, solver = solver)
  invisible(obj)
} 

setMethod("initialize", signature(.Object = "polyPdmpModel"),
          function(.Object, 
                   discDomain, dynpolys, ratepolys, jumpfunc, parms, init, 
                   dynfunc = NULL, ratefunc = NULL, out = NULL, label = character(0), descr = character(0),
                   dynsprays = NULL, ratesprays = NULL,
                   times = c(from = 0, to = 10, by = 1), 
                   solver = "lsodar", ...){
            
            if(!is.null(dynfunc)) warning("The input for 'dynfunc' is ignored, because dynfunc is automatically specified by 'dynpolys'.")
            if(!is.null(ratefunc)) warning("The input for 'ratefunc' is ignored, because ratefunc is automatically specified by 'ratepolys'.")
            if(!is.null(ratesprays)) warning("The input for 'ratesprays' is ignored, because ratesprays is automatically specified by 'ratepolys' and 'parms'.")
            if(!is.null(dynsprays)) warning("The input for 'dynsprays' is ignored, because dynsprays is automatically specified by 'dynpolys' and 'parms'.")
            
            .Object@discDomain <- discDomain # needed vor definition of dynpolys and ratepolys
            .Object@init <- init             # needed for definition of dynpolys and ratepolys
            .Object@parms <- parms           # needed for definition of dynpolys and ratepolys
            
            dynpolys(.Object) <- dynpolys    # also specifies the slots "dynfunc" and ".evaluatedDynpolys"
            ratepolys(.Object) <- ratepolys  # also specifies the slots "ratefunc" and ".evaluatedRatepolys"
            
            validObject(.Object)
            callNextMethod(.Object, jumpfunc = jumpfunc,  times = times, solver = solver, out = out, label = label, descr = descr, ...)
            }) 

setValidity("polyPdmpModel", function(object){
  return(TRUE)
})
