#======== todo =================================================================

setClass("polyPdmpModel", 
         slots = list(dynpolys = "call", ratepolys = "call", 
                      dynsprays = "list", ratesprays = "list"), 
         contains = "pdmpModel")

polyPdmpModel <- function(obj = NULL, .cache = new.env(),
                          discStates, dynpolys, ratepolys, jumpfunc, 
                          parms, init, dynfunc = NULL, ratefunc = NULL, 
                          out = NULL, descr = character(0),
                          dynsprays = NULL, ratesprays = NULL,
                          times = c(from = 0, to = 10, by = 1), 
                          solver = "lsodar", ...){
  # see function "initialize" for detailed code
  obj <- new(Class = "polyPdmpModel", discStates = discStates, 
             dynpolys = dynpolys, dynfunc = dynfunc, jumpfunc = jumpfunc, 
             ratepolys = ratepolys, ratefunc = ratefunc, dynsprays = dynsprays, 
             ratesprays = ratesprays, descr = descr, solver = solver,
             parms = parms, out = out, init = init, times = times)
  invisible(obj)
} 

setMethod("initialize", signature(.Object = "polyPdmpModel"),
          function(.Object, 
                   discStates, dynpolys, ratepolys, jumpfunc, parms, init, 
                   dynfunc = NULL, ratefunc = NULL, out = NULL, 
                   descr = character(0), dynsprays = NULL, ratesprays = NULL,
                   times = c(from = 0, to = 10, by = 1), 
                   solver = "lsodar", ...){
            
            if(!is.null(dynfunc)) 
              warning("The input for 'dynfunc' is ignored, because dynfunc is 
                      automatically specified by 'dynpolys'.")
            if(!is.null(ratefunc)) 
              warning("The input for 'ratefunc' is ignored, because ratefunc 
                      is automatically specified by 'ratepolys'.")
            if(!is.null(ratesprays)) 
              warning("The input for 'ratesprays' is ignored, because 
                      ratesprays is automatically specified by 'ratepolys' 
                      and 'parms'.")
            if(!is.null(dynsprays)) 
              warning("The input for 'dynsprays' is gnored, because dynsprays 
                      is automatically specified by 'dynpolys' and 'parms'.")
            
            # needed for definition of dynpolys and ratepolys:
            .Object@discStates <- discStates 
            .Object@init <- init             
            .Object@parms <- parms           
            
            dynpolys(.Object) <- dynpolys    # also specifies the slots "dynfunc" and ".evaluatedDynpolys"
            ratepolys(.Object) <- ratepolys  # also specifies the slots "ratefunc" and ".evaluatedRatepolys"
            
            validObject(.Object)
            callNextMethod(.Object, jumpfunc = jumpfunc,  times = times, 
                           solver = solver, out = out, descr = descr, ...)
            }) 

setValidity("polyPdmpModel", function(object){
  return(TRUE)
})
