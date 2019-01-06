#======== todo =================================================================
#t2 hier Dokumentation nicht doppelt. Geht das auch?
#t3 validity schreiben (was kommt da genau rein?, siehe validObject(.Object))
#t2 arity(spray) == length(obj@init) für dynpolys und ratepolys
#s3 Bei polyPdmpModel<- wird .cache = new.env() gesetzt, bei pdmpModel<- nicht. 
#   Warum? Im Internet habe ich bisher nicht viel dazu gefunden. Vllt wegen initialize?

#' Class polyPdmpModel
#' 
#' An S4 class to represent polynomial piecewisc deterministic markov processes
#' (polynomial PDMPs). These processes are PDMPs with polynomial rate functions
#' and polynomial dynamics. This makes it possible to approximate the moments of
#' the process without the need of simulation, see \code{\link{momApp}}. On the
#' downside, simulations take much longer than simulations of the same PDMP 
#' represented as \code{pdmpModel}. \cr  
#' The class is based on the \code{\link[pdmpsim]{pdmpModel}} class of package
#' \pkg{pdmpsim} but introduces additional slots that replace the slots
#' \code{dynfunc} and \code{ratefunc}, namely \code{dynpolys}, \code{dynsprays},
#' \code{ratepolys} and \code{ratesprays}. Only slots \code{dynpolys} and
#' \code{ratepolys} have to be defined by the user. The other slots (including
#' the still existing \code{dynfunc} and \code{ratefunc}) will be set
#' automatically. To represent polynomials in \code{R}, package \pkg{spray} is
#' used, which stores the coefficient matrices of polynomials as sparse arrays.
#' 
#' @slot dynpolys an object of class \code{language}, more precisely a quoted
#'   list. Every element of this list contains the ODE for one of the continous
#'   variables. These ODEs are given as \code{spray} objects and can be defined
#'   in two different ways. See section dynpolys for further information.
#' @slot ratepolys an object of class \code{language}, more precisely a quoted
#'   list.This list contains the rates for every jumptype. Every element is
#'   itself a list, its length is determined by the number of different values
#'   of the discrete variable (\code{length(discStates(obj)[[1]])}). See section
#'   ratepolys for further information.
#' @slot dynsprays a nested list of spray objects. This slot is generated
#'   automatically out of slot \code{dynpolys}. Do not edit it by hand.
#' @slot ratesprays a nested list of spray objects. This slot is generated
#'   automatically out of slot \code{ratepolys}. Do not edit it by hand.
#' @note Currently, all methods only work for PDMPs with a single discrete
#'   variable, but it is always possible to redefine a PDMP with multiple
#'   discrete variables in a way that only one discrete variable is necessary.
#' 
#' @section Ratepolys:
#' Slot ratepolys is a quted list. The length of the list determines the number
#' of existing jumptypes (jumptypes are used in slot \code{jumpfunc} to
#' determine the next state the process jumps too). It contains the rates
#' determining the probability of a jump to be of type 1, 2, and so on. The
#' rates are polynomials that usually depend on the value of the discrete
#' variable at the time when the jump occurs. They are given as a list of spray
#' objects or numbers, for every possible discrete state separately (in the same
#' order as the states are given in slot \code{discStates}). \cr
#' 
#' For example, let's assume that we have three different jumptypes and a
#' discrete variable \code{d} that can take the values 0 or 1 (so slot
#' \code{discStates} would be defined as \code{list(d = 0:1)}). Then slot
#' \code{ratepolys} will be given as follows:
#' \preformatted{quote(list(
#'   list(rate of jtype 1 with d = 0, rate of jtype 1 with d = 1),
#'   list(rate of jtype 2 with d = 0, rate of jtype 2 with d = 1),
#'   list(rate of jtype 3 with d = 0, rate of jtype 3 with d = 1),
#' ))}
#' Every \code{rate of jtype i with d = j} is a spray object, for example 
#' \code{linear(1:3)}.
#' 
#' @section Dynpolys:
#' Slot dynpolys is a quoted list. The length of the list equals the number of
#' continuous variables of the model. Every element contains the ODE for one of
#' the continuous variables, given in the same order as in slot \code{init} and
#' given as spray objects. \cr The ODEs usually depend on the discrete variable,
#' lets call it \code{d}. An example would be \eqn{\frac{df}{dt} = -3f}{df/dt =
#' -3f} if \eqn{d = 0} and  \eqn{\frac{df}{dt} = -3f + 1}{df/dt = -3f + 1} if
#' \eqn{d = 1}. Here, \code{f} is the only continous variable of the model.
#' Every ODE can be defined in two different ways:
#' 
#' \itemize{
#' \item variant 1: \cr
#' a list whose length equals the number of different states for d. Every entry is a 
#' polynomial with d taking a fixed value. The order of the entries corresponds
#' to the order of discrete states given in slot \code{discStates}.
#' \item variant 2: \cr
#' a list of two variables:
#' \itemize{
#' \item a variable \code{overall} which contains a polynomial that is independent 
#' of fixed values for \code{d} and is therefore the same for all states 
#' (attention: \code{overall} can contain \code{d} as a formal variable),
#' \item a variable \code{specific} which contains the rest of the term.
#' This is a list of the same form as in variant 1.
#' }}
#' In our example, we get
#' \itemize{
#' \item Variant 1: \cr
#' \code{quote(list(list(-3*lone(1,2), -3*lone(1,2) + 1)))}
#' \item Variant 2: \cr
#' \code{quote(list(list(overall = linear(c(-3, 0)), specific = list(0, 1))))}
#' \item Variant 2 with one formula: \cr
#' \code{quote(list(list(overall = linear(c(-3, 1))))}
#' }
#' 
#' The last variant is only possible, because in this example we have the
#' possiblity to write both ODEs in one formula: \eqn{\frac{df}{dt} = -3f +
#' d}{df/dt = -3f + d}.
#' @example /inst/examples/polyPdmpModel.R
#' @aliases polyPdmpModel polypdmpmodel polypdmp polyPdmp
#' @importFrom pdmpsim pdmpModel
#' @importFrom methods new
#' @export
setClass("polyPdmpModel", 
         slots = list(dynpolys = "call", ratepolys = "call", 
                      dynsprays = "list", ratesprays = "list"), 
         contains = "pdmpModel")

polyPdmpModel <- function(obj = NULL, .cache = new.env(),
                          discStates, dynpolys, ratepolys, jumpfunc, 
                          parms = c(0), init, dynfunc = NULL, ratefunc = NULL, 
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

#' @importFrom methods validObject callNextMethod
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
