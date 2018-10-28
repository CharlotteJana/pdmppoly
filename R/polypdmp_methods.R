#======== todo =================================================================
#t3 länge der liste mit anzahl von jtypes vergleichen (if(is.null(where)) ...)

#' @include polypdmp_class.R spray_extension.R
NULL

#' Redefine Slot \code{ratepolys}
#'
#' This function is used internally for function \code{ratepolys<-} that 
#' defines slot \code{ratepolys} and during initialization of a
#' \code{\link{polyPdmpModel}} object. It replaces parameters stored in slot
#' \code{parms} by its values and coerces all numbers to \code{spray} objects.
#' @param x an object of type language as described in section Ratepolys in
#'   \code{\link{polyPdmpModel}}.
#' @param obj an object of class polyPdmpModel.
#' @param where this parameter is necessary because \code{redefineRatepolys} 
#' is a recursive function. Changing it will have no effect on the result.
#' @return an modified object of type language as described in section 
#' Ratepolys in \code{\link{polyPdmpModel}}.
#' @examples 
#' data(simplePoly)
#' newRates <- quote(list(list(1, 2, 3), 
#'                   list(lone(1,2), linear(3:4), product(1:2))))
#' redefineRatepolys(newRates, simplePoly)                  
#' @importFrom spray is.spray
#' @export
redefineRatepolys <- function(x, obj, where = NULL){
  
  dim <- length(obj@init)
  eva <- with(as.list(obj@parms), eval(x))

  if(is.numeric(eva)) 
    bquote(.(x)*one(.(dim))) # turn numbers into spray-objects
  else if(is.spray(eva))
    x
  else if(is.list(eva)) # call redefineRatepolys for all list entries
    as.call(lapply(x, redefineRatepolys, obj = obj, where = parent.frame()))
  else if(is.name(x) && identical(x, as.name("list"))) 
    x
  else 
    stop("Input '", x, "' in 'ratepolys' is not correct.")
}

#' Redefine Slot \code{dynpolys}
#'
#' This function is used internally for function \code{dynpolys<-} that defines
#' slot \code{dynpolys} and during initialization of a
#' \code{\link{polyPdmpModel}} object. It replaces parameters stored in slot
#' \code{parms} by its values and coerces all numbers to \code{spray} objects.
#' It coerces all list elements that are given in variant 2 to spray objects
#' defined as in variant 1 (see section Dynpolys in \code{\link{polyPdmpModel}}
#' for a description of the two variants).
#' @param x an object of type language as described in section Dynpolys in
#'   \code{\link{polyPdmpModel}}.
#' @param obj an object of class polyPdmpModel.
#' @param where this parameter is necessary because \code{redefineDynpolys} 
#' is a recursive function. Changing it will have no effect on the result.
#' @param overall this parameter is necessary because \code{redefineDynpolys} 
#' is a recursive function. Please do not change its value!
#' @return a modified object of type language as described in section 
#' Dynpolys in \code{\link{polyPdmpModel}}.
#' @note This function only works for one discrete variable
#' @importFrom spray is.zero is.spray
#' @export
redefineDynpolys <- function(x, obj, where = NULL, overall = 0){
  
  eva <- with(as.list(obj@parms), eval(x))
  eoa <- with(as.list(obj@parms), eval(overall))
  dim <- length(obj@init)
  
  # if 'overall' == number: turn it into a 'spray'-object
  if(is.numeric(eoa)) { 
    overall <- bquote(.(overall)*one(.(dim)));
    eoa <- with(as.list(obj@parms), eval(overall))
  }
  
  # if x == value (number or 'spray'): 
  # add 'overall' to x and turn x into a 'spray'-object (when necessary)
  if(is.numeric(eva)) {
    if(is.zero(eoa)) bquote(.(x)*one(.(dim))) 
    else  bquote(.(overall)+.(x)*one(.(dim)))
  }
  else if(is.spray(eva)){
    if(is.zero(eoa)) x
    else  bquote(.(overall)+.(x))
  }
  
  # if x == 'list':  rerun 'redefineDynpolys'
  else if(is.name(x) && identical(x, as.name("list"))) 
    x
  else if(is.list(eva)) {
    if(!is.null(x$overall) & is.null(x$specific)) {
      overall <- x$overall
      zeroList <- as.list(rep(0, length(obj@discStates[[1]])+1));
      zeroList[[1]] <- "list"
      x <- do.call(call, zeroList)
    } 
    if(!is.null(x$overall)){
      overall <- x$overall
      x$overall <- NULL
    }
    if(!is.null(x$specific))
      x <- x$specific
    
    as.call(lapply(x, redefineDynpolys, obj, where = parent.frame(), overall))
  }
  
  # if x == something else: give error
  else {stop("Input ", x, " for 'dynpolys' is not correct.")}
}  

#' get Index of variable 'var' in vector 'vect'
#' 
#' This method is used internally. If no index is found with \code{\link{which}}
#' (this can happen because of numerical issues), the less efficient method
#' \code{\link{all.equal}} is used to check for "nearly" equality.
#' @param vect vector
#' @param var element that one expects to be part of \code{vect}
#' @examples
#' getIndex("c", letters)
#' \donttest{getIndex(5.0000001, 2:6)} # error
#' getIndex(5.00000001, 2:6) # difference is small enough to be found by all.equal
#' @export
getIndex <- function(var, vect){
  index <- which(vect == var) 
  if(length(index) == 0){
    #print("all.equal was used")
    comp <- mapply(function(x) {isTRUE(all.equal(x, var, check.names = FALSE))}, vect) 
    index <- which(comp == TRUE)
  }
  if(length(index) == 0)
    stop("Element ", var, " not present in vector ", paste(vect, collapse = ", "))
  return(index)
}

#' Convert ratesprays into a matrix
#' 
#' This method computes a 'ratematrix' out of slot \code{ratesprays}.
#' The i,j - entry will then be accessible via \code{ratematrix[[i]][[j]]}
#' and contain the rate from discrete state \code{i} to discrete state \code{j}.
#' Note that the returned object is not of class \code{matrix}, but is a
#' nested list of \code{spray} objects.
#' @param obj object of class \code{\link{polyPdmpModel}}.
#' @note this method only works for one discrete variable
#' @return a two dimensional list of \code{spray} objects.
#' @export
ratespraysToMatrix <- function(obj){
  
  n <- length(obj@init) - length(obj@discStates) # number of continuous variables
  nj <- length(obj@ratesprays) # number of jumptypes
  k <- length(obj@discStates[[1]]) # number of different discrete states
  z <- obj@init
  discIndex <- which(names(obj@init) == names(obj@discStates[1]))
  
  ratematrix <- rep(list(rep(list(NULL), k)), k)
  for(jtype in 1:nj){ 
    for(i in 1:k){
      oldDiscVar <- obj@discStates[[1]][[i]]
      z[discIndex] <- oldDiscVar
      for(j in 1:k){
        newDiscVar <- obj@discStates[[1]][[j]]
        if(obj@jumpfunc(1, z, obj@parms, jtype)[discIndex] == newDiscVar){
          ratematrix[[i]][[j]] <- obj@ratesprays[[jtype]][[i]]
        }
      } 
    }
  }
  return(ratematrix)
}

#' Change discrete variable in a polynomial
#' 
#' This function takes a polynomial (represented by a \code{spray} object),
#' uses polyPdmpModel \code{obj} to identify the discrete variable and
#' substitutes this variable in the polynomial by several indicator
#' variables. The number of variables depends on the different possible
#' states that the discrete variable can take. The method is used internally
#' for computing the \code{\link{generator}} of a \code{polyPdmpModel}.
#' 
#' @details
#' This function is best explained in an example. Lets assume that
#' we have one continous variable \code{f} and one discrete variable \code{d}
#' and that they are given in order \code{f, d} in variable \code{init(obj)}.
#' Let's have \code{discStates(obj) = list(d = -1:1)}, so \code{d} can take
#' the values -1, 0 or 1, respectively. A polynomial of arity 2 would then be
#' transformed into a polynomial of arity 4, where the second variable (that 
#' stands for \code{d}) is replaced by three indicator variables 
#' \code{d1, d2, d3}, that can only take values 0 or 1, depending on the
#' value of the original \code{d}. If we take for example the polynomial
#'  \eqn{5f^2d^3 + f}{5f²d³ + f}, method \code{blowupSpray} changes it into
#' \deqn{
#' 5f^2\cdot(-1)^3\cdot d1 + 5f^2\cdot 0^3 \cdot d2 + 5f^2\cdot 1^3\cdot d3 + f 
#' = -5f^2d1 + 5f^2d3 + f.}{
#' 5f²⋅(-1)³⋅d1 + 5f²⋅0³⋅d2 + 5f²⋅1³⋅d3 + f = -5f²d1 + 5f²d3 + f.}
#' If the order of \code{f} and \code{d} are reversed (i.e. \code{init(obj) =
#' c(d = ..., f = ...)}), the additional variables \code{d1, d2, d3} would 
#' still be added at the end, leaving \code{f} to be the first variable.
#' 
#' @examples
#' # the example discussed in section details:
#' library(spray)
#' data("simplePoly") 
#' blowupSpray(simplePoly, 5*product(2:3) + lone(1, 2))
#' 
#' # with other order in slot init:
#' init(simplePoly) <- c(d = 1, f = 2)
#' blowupSpray(simplePoly, 5*product(3:2) + lone(2, 1))  
#' 
#' @param obj object of class \code{\link{polyPdmpModel}}.
#' @param spray object of class \code{spray}. This is the 
#' polynomial that shall be modified.
#' @note This method only works for one discrete variable
#' @importFrom spray is.zero subs lone
#' @export
blowupSpray <- function(obj, spray){
  
  if(is.zero(spray)) return(spray)
  
  n <- length(obj@init) # number of continous variables
  k <- length(obj@discStates[[1]]) # number of different discrete states
  
  stopifnot(k > 0)
  stopifnot(arity(spray) == n)
  
  #spray <- increase_arity(spray, n + 1:(k-1))
  
  # part of spray that is independent of discrete variable:
  discIndex <- which(names(obj@init) == names(obj@discStates[1]))
  c <- subs(spray, discIndex, 0)
  if(is.zero(c)){
    coreShort <- 0*lone(1, n)
    core <- 0*lone(1, n-1+k)
  }
  else{
    coreShort <- increase_arity(c, discIndex)
    core <- increase_arity(c, n + 0:(k-1))
  }
  
  # add the specific parts (depend on discrete variable):
  blowedSpray <- core
  for(i in 1:k){
    discVar <- obj@discStates[[1]][i]
    if(!is.zero(spray - coreShort)){
      specific <- subs(spray - coreShort, discIndex, discVar)
      if(!is.zero(specific)){
        specific <- increase_arity(specific,  n + 0:(k-1))*lone(n-1+i, n-1+k)
        blowedSpray <- blowedSpray + specific
      }
    }
  }
  return(blowedSpray)
}


##### output methods ####

#' @importFrom methods slot slotNames
#' @export
setMethod(f = "print",
          signature = "polyPdmpModel",
          definition = function(x, all = FALSE, ...){
              if(all){
                print.default(x, all = TRUE, ...)
              }
              else {
                cat("An S4-object of class", class(x)[1], "\n\n")
                slotnames <- slotNames(x)
                for (slotname in slotnames) {
                  slotcontent <- slot(x, slotname)
                  if (!is.null(slotcontent)) {
                    if (slotname == "main" || 
                        slotname == "dynfunc" ||
                        slotname == "dynsprays" ||
                        slotname == "ratesprays" ||
                        slotname == "ratefunc") next
                    cat("Slot ", dQuote(slotname), ":\n", sep = "")
                    if (slotname == "out")    cat("  outputs exist ...\n")
                    else if (slotname == "parms")  str(slotcontent)
                    else {
                      print(slotcontent)
                    }
                    cat("\n")
                  }
                }
                cat("Hint: use print(x, all=TRUE) to see all polyPdmpModel slots.\n")
              }
          }
)
