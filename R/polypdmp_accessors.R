setGeneric("dynpolys", function(obj, ...) standardGeneric("dynpolys"))
setGeneric("dynpolys<-", function(obj, value) standardGeneric("dynpolys<-"))
setGeneric("ratepolys", function(obj, ...) standardGeneric("ratepolys"))
setGeneric("ratepolys<-", function(obj, value) standardGeneric("ratepolys<-"))
setGeneric("ratesprays", function(obj, ...) standardGeneric("ratesprays"))
setGeneric("ratesprays<-", function(obj, value) standardGeneric("ratesprays<-"))
setGeneric("dynsprays", function(obj, ...) standardGeneric("dynsprays"))
setGeneric("dynsprays<-", function(obj, value) standardGeneric("dynsprays<-"))

setMethod("dynpolys", "polyPdmpModel", function(obj, ...) obj@dynpolys)
setMethod("dynpolys<-", "polyPdmpModel", function(obj, value){
  
  ### dynpolys is a list with length = (number of continuous variables)
  ### every entry contains the ode for one of the continuous variables
  ### these terms depend on the discrete parameter θ and can be defined in two different variants:
  ### variant 1: a list of length = (number of different states for θ), every entry is a polynomial with θ fixed
  ### variant 2: a list of two variables:
  ###            a variable 'overall' which contains a polynomial that is independent of fixed values for θ 
  ###              and is therefore the same for all states (attention: 'overall' can contain θ as a formal variable),
  ###            a variable 'specific' which contains the rest of the term, it is a list of the same form as in variant 1

  obj@dynpolys <- value
  
  redefined <- redefineDynpolys(value, obj)
  # check if 'value' is defined correctly,
  # all entrys of 'value' are coerced to the form 'variant 1',
  # numbers are coerced to 'spray'-objects.
  # (This is still an object of type 'language')
  
  obj@dynsprays <- with(as.list(obj@parms), eval(redefined))
  # evaluate dynpolys, the values of the parameters are taken from obj@parms,
  # (This is a nested list of sprays)
  
  obj@dynfunc <- function(t, z, parms = obj@parms){
    discDomainIndex <- getIndex(z[length(z)], obj@discDomain) # index of discDomain that corresponds to the current value of the discrete variable
    if(!identical(parms,obj@parms)) {
      stop("please redefine the slot 'parms' (by using 'parms<-).")
    }
    else dynsprays <- obj@dynsprays
    
    funcs <- lapply(dynsprays, function(x) lapply(x, as.function.spray))
    funcs <- lapply(funcs, function(list) list[[discDomainIndex]]) # pick the right sprays out of dynpolys (one for every continous variable)
    dz <- sapply(funcs, function(f) unname(do.call(f, list(z)))) # apply these functions to z
    return(c(dz,0))
  }
  # turn all 'dynsprays' into functions(z),
  # pick the right ones (depending on discvar = z[n+1]) and apply them to z.
  # (This is a vector of values with length = number of continuous variables).
  # Note: there is no real dependence of parms, because parms and obj@parms have to be equal.

  obj@out <- NULL
  invisible(obj)
})
setMethod("ratepolys", "polyPdmpModel", function(obj, ...) obj@ratepolys)
setMethod("ratepolys<-", "polyPdmpModel", function(obj, value){
  
  obj@ratepolys <- value
    
  redefined <- redefineRatepolys(value, obj)
  # check if 'value' is defined correctly,
  # numbers are coerced to 'spray'-objects.
  # (This is still an object of type 'language')
  
  obj@ratesprays <-  with(as.list(obj@parms), eval(redefined))
  # evaluate ratepolys, the values of the parameters are taken from obj@parms,
  # (This is a nested list of sprays)
  
  obj@ratefunc <- function(t, z, parms = obj@parms){
    discDomainIndex <- getIndex(z[length(z)], obj@discDomain) # index of discDomain that corresponds to the current value of the discrete variable
    if(!identical(parms,obj@parms)){ 
      stop("please redefine the slot 'parms' (by using 'parms<-).")
    }                     
    else ratesprays <- obj@ratesprays
    
    funcs <-lapply(ratesprays, function(x) lapply(x, as.function.spray))
    funcs <- lapply(funcs, function(list) list[[discDomainIndex]]) # pick the right sprays out of ratepolys (one for every jumptype)
    return(sapply(funcs, function(f) unname(do.call(f, list(z))))) # apply these functions to z
  }
  # turn all 'ratesprays' into functions(z),
  # pick the right ones (depending on discvar = z[n+1]) and apply them to z.
  # (This is a vector of values with length = number of jumptypes).
  # Note: there is no real dependence of parms, because parms and obj@parms have to be equal.
  
  obj@out <- NULL
  invisible(obj)
})
setMethod("parms<-", "polyPdmpModel", function(obj, value){
  obj@parms <- value
  obj@dynsprays <- with(as.list(value), eval(redefineDynpolys(obj@dynpolys, obj)))
  obj@ratesprays <- with(as.list(value), eval(redefineRatepolys(obj@ratepolys, obj)))
  out(obj) <- NULL
  invisible(obj)
})

setMethod("dynfunc", "polyPdmpModel", function(obj, ...) obj@dynfunc)
setMethod("dynfunc<-", "polyPdmpModel", function(obj, value) {
  stop("The slot 'dynfunc' is automatically generated using 'dynpolys'. Please modify 'dynpolys' instead.")
})
setMethod("ratefunc", "polyPdmpModel", function(obj, ...) obj@ratefunc)
setMethod("ratefunc<-", "polyPdmpModel", function(obj, value) {
  stop("The slot 'ratefunc' is automatically generated using 'rateploys'. Please modify 'ratepolys' instead.")
})
setMethod("ratesprays", "polyPdmpModel", function(obj, ...) obj@ratesprays)
setMethod("ratesprays<-", "polyPdmpModel", function(obj, value) {
  stop("The slot 'ratesprays' is automatically generated using 'ratepolys' and 'parms'. Please modify these slots instead.")
})
setMethod("dynsprays", "polyPdmpModel", function(obj, ...) obj@dynsprays)
setMethod("dynsprays<-", "polyPdmpModel", function(obj, value) {
  stop("The slot 'dynsprays' is automatically generated using 'dynpolys' and 'parms'. Please modify these slots instead.")
})
