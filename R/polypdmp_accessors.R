#======== todo =================================================================

#' Accessor functions for Class polyPdmpModel
#' 
#' Get or set the slot values for an object of class
#' \code{\link{polyPdmpModel}}. There are no methods for setting new values of
#' the slots \code{dynsprays}, \code{dynfunc}, \code{ratesprays} and 
#' \code{ratefunc}. Their value can be changed by setting a new value to
#' \code{dynpolys}, \code{ratepolys} or \code{parms}. Only these slots
#' should be accessed by the user. Slot \code{out} will be set to NULL, if
#' new values for the other slots are provided.
#' 
#' It is also possible to modify the slots directly,
#' e.g. slot \code{dynsprays} of a model \code{model} with \code{model@dynsprays}, 
#' but this is not recommended! The slots have to fit together and there is
#' no method warning you of possible misleading results.
#' 
#' @param obj an object of class \code{\link{polyPdmpModel}}
#' @param value the value that shall be set
#' 
#' @examples
#' data("genePolyK")
#' parms(genePolyK)["b"] <- 3 
#' dynsprays(genePolyK) # different value than before
#' 
#' @include polypdmp_class.R polypdmp_methods.R
#' @importMethodsFrom pdmpsim "dynfunc<-" "ratefunc<-" "parms<-"
#' @name polypdmp-accessors
NULL

#========= setGenerics ===========

#' @rdname polypdmp-accessors
#' @export
setGeneric("dynpolys", function(obj) standardGeneric("dynpolys"))
#' @rdname polypdmp-accessors
#' @export
setGeneric("dynpolys<-", function(obj, value) standardGeneric("dynpolys<-"))
#' @rdname polypdmp-accessors
#' @export
setGeneric("ratepolys", function(obj) standardGeneric("ratepolys"))
#' @rdname polypdmp-accessors
#' @export
setGeneric("ratepolys<-", function(obj, value) standardGeneric("ratepolys<-"))
#' @rdname polypdmp-accessors
#' @export
setGeneric("ratesprays", function(obj) standardGeneric("ratesprays"))
#' @rdname polypdmp-accessors
#' @export
setGeneric("ratesprays<-", function(obj, value) standardGeneric("ratesprays<-"))
#' @rdname polypdmp-accessors
#' @export
setGeneric("dynsprays", function(obj) standardGeneric("dynsprays"))
#' @rdname polypdmp-accessors
#' @export
setGeneric("dynsprays<-", function(obj, value) standardGeneric("dynsprays<-"))
# #' @rdname pdmp-accessors
# #' @export
# setGeneric("parms<-", function(obj, value) standardGeneric("parms<-"))

#======= Getters ==========

#' @rdname polypdmp-accessors
#' @export
setMethod("dynpolys", "polyPdmpModel", function(obj) obj@dynpolys)
#' @rdname polypdmp-accessors
#' @export
setMethod("ratepolys", "polyPdmpModel", function(obj) obj@ratepolys)
#' @rdname polypdmp-accessors
#' @export
setMethod("ratesprays", "polyPdmpModel", function(obj) obj@ratesprays)
#' @rdname polypdmp-accessors
#' @export
setMethod("dynsprays", "polyPdmpModel", function(obj) obj@dynsprays)

#======= Setters ==========

#' @note only works for one discrete variable
#' @rdname polypdmp-accessors
#' @importFrom spray as.function.spray
#' @export
setMethod("dynpolys<-", "polyPdmpModel", function(obj, value){
  
  obj@dynpolys <- value
  
  redefined <- redefineDynpolys(value, obj)
  # check if 'value' is defined correctly,
  # all entrys of 'value' are coerced to the form 'variant 1',
  # numbers are coerced to 'spray'-objects.
  # (This is still an object of type 'language')
  
  obj@dynsprays <- with(as.list(obj@parms), eval(redefined))
  # evaluate dynpolys, the values of the parameters are taken from obj@parms,
  # (This is a nested list of sprays)
  
  obj@dynfunc <- function(t, x, parms = obj@parms){
    discName <- names(obj@discStates)
    discStatesIndex <- getIndex(x[discName], obj@discStates[[1]]) # index of discStates that corresponds to the current value of the discrete variable
    if(!identical(parms, obj@parms)) {
      stop("please redefine the slot 'parms'.")
    }
    else dynsprays <- obj@dynsprays
    
    funcs <- lapply(dynsprays, function(y) lapply(y, as.function.spray))
    funcs <- lapply(funcs, function(list) list[[discStatesIndex]]) # pick the right sprays out of dynpolys (one for every continous variable)
    dx <- sapply(funcs, function(f) unname(do.call(f, list(x)))) # apply these functions to x
    contNames <- names(obj@init)[!names(obj@init) %in% discName]
    dynamics <- rep(0, length(obj@init)) %>% stats::setNames(names(obj@init))
    dynamics[contNames] <- dx
    return(dynamics)
  }
  # turn all 'dynsprays' into functions(x),
  # pick the right ones (depending on discvar) and apply them to x.
  # (This is a vector of values with length = number of continuous variables).
  # Note: there is no real dependence of parms, because parms and obj@parms have to be equal.

  obj@out <- NULL
  invisible(obj)
})

#' @note only works for one discrete variable
#' @rdname polypdmp-accessors
#' @importFrom spray as.function.spray
#' @export
setMethod("ratepolys<-", "polyPdmpModel", function(obj, value){
  
  obj@ratepolys <- value
    
  redefined <- redefineRatepolys(value, obj)
  # check if 'value' is defined correctly,
  # numbers are coerced to 'spray'-objects.
  # (This is still an object of type 'language')
  
  obj@ratesprays <-  with(as.list(obj@parms), eval(redefined))
  # evaluate ratepolys, the values of the parameters are taken from obj@parms,
  # (This is a nested list of sprays)
  
  obj@ratefunc <- function(t, x, parms = obj@parms){
    discStatesIndex <- getIndex(x[names(obj@discStates)], obj@discStates[[1]]) # index of discStates that corresponds to the current value of the discrete variable
    if(!identical(parms,obj@parms)){ 
      stop("please redefine the slot 'parms' (by using 'parms<-).")
    }                     
    else ratesprays <- obj@ratesprays
    
    funcs <-lapply(ratesprays, function(y) lapply(y, as.function.spray))
    funcs <- lapply(funcs, function(list) list[[discStatesIndex]]) # pick the right sprays out of ratepolys (one for every jumptype)
    return(sapply(funcs, function(f) unname(do.call(f, list(x))))) # apply these functions to x
  }
  # turn all 'ratesprays' into functions(x),
  # pick the right ones (depending on discvar = x[n+1]) and apply them to x.
  # (This is a vector of values with length = number of jumptypes).
  # Note: there is no real dependence of parms, because parms and obj@parms have to be equal.
  
  obj@out <- NULL
  invisible(obj)
})

#' @rdname polypdmp-accessors
#' @export
setMethod("parms<-", "polyPdmpModel", function(obj, value){
  obj@parms <- value
  obj@dynsprays <- with(as.list(value), eval(redefineDynpolys(obj@dynpolys, obj)))
  obj@ratesprays <- with(as.list(value), eval(redefineRatepolys(obj@ratepolys, obj)))
  obj@out <- NULL
  invisible(obj)
})

#======= Setters that throw errors ==========

#' @rdname polypdmp-accessors
#' @export
setMethod("dynfunc<-", "polyPdmpModel", function(obj, value) {
  stop("The slot 'dynfunc' is automatically generated using 'dynpolys'. 
       Please modify 'dynpolys' instead.")
})
#' @rdname polypdmp-accessors
#' @export
setMethod("ratefunc<-", "polyPdmpModel", function(obj, value) {
  stop("The slot 'ratefunc' is automatically generated using 'rateploys'. 
       Please modify 'ratepolys' instead.")
})
#' @rdname polypdmp-accessors
#' @export
setMethod("ratesprays<-", "polyPdmpModel", function(obj, value) {
  stop("The slot 'ratesprays' is automatically generated using 'ratepolys' 
       and 'parms'. Please modify these slots instead.")
})
#' @rdname polypdmp-accessors
#' @export
setMethod("dynsprays<-", "polyPdmpModel", function(obj, value) {
  stop("The slot 'dynsprays' is automatically generated using 'dynpolys' 
       and 'parms'. Please modify these slots instead.")
})
