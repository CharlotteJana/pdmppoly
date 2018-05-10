#' @include polypdmp_class.R
NULL

# needed for "ratepolys←" and "dynpolys←" :
#' @importFrom spray is.spray
redefineRatepolys <- function(x, obj, where = parent.frame()){
  # check if 'x' is defined correctly, numbers are coerced to 'spray'-objects.
  # (This is still an object of type 'language')
  
  eva <- with(as.list(obj@parms), eval(x)) #  only the names of the parameters are important
  dim <- length(obj@init)
  if(is.numeric(eva)) {bquote(.(x)*one(.(dim)))} # turn numbers into 'spray'-objects
  else if(is.spray(eva)) {x}
  else if(is.list(eva)) {as.call(lapply(x, redefineRatepolys, obj = obj, where = where))}
  else if(is.name(x) && identical(x, as.name("list"))) {x}
  else { stop("Input '", x, "' in 'ratepolys' is not correct.")}
}

#' @note only works for one discrete variable
#' @importFrom spray is.zero is.spray
redefineDynpolys <- function(x, obj, where = parent.frame(), overall = 0){
  # check if 'x' is defined correctly, numbers are coerced to 'spray'-objects,
  # all entrys of 'x' are coerced to the form 'variant 1' (description: see "dynpolys←")
  # (This is still an object of type 'language')
  
  eva <- with(as.list(obj@parms), eval(x)) #  only the names of the parameters are important
  eoa <- with(as.list(obj@parms), eval(overall))
  dim <- length(obj@init)
  
  # if 'overall' == number: turn it into a 'spray'-object
  if(is.numeric(eoa)) { 
    overall <- bquote(.(overall)*one(.(dim)));
    eoa <- with(as.list(obj@parms), eval(overall))
  }
  
  # if x == value (number or 'spray'): add 'overall' to x and turn x into a 'spray'-object (when necessary)
  if(is.numeric(eva)) {
    if(is.zero(eoa)) bquote(.(x)*one(.(dim))) # if 'overall' == 0: don't add anything, turn x into 'spray'
    else  bquote(.(overall)+.(x)*one(.(dim)))
  }
  else if(is.spray(eva)){
    if(is.zero(eoa)) x
    else  bquote(.(overall)+.(x))
  }
  
  # if x == 'list':  rerun 'redefineDyn'
  else if(is.name(x) && identical(x, as.name("list"))) {x}
  else if(is.list(eva)) {
    if(!is.null(x$overall) & is.null(x$specific)) {
      overall <- x$overall
      zeroList <- as.list(rep(0, length(obj@discStates[[1]])+1));
      zeroList[[1]] <- "list"
      x <- do.call(call, zeroList)
    } 
    if(!is.null(x$overall)) {overall <- x$overall; x$overall <- NULL}
    if(!is.null(x$specific)) {x <- x$specific}
    as.call(lapply(x, redefineDynpolys, obj = obj, where = where, overall = overall))
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
#' getIndex(5.0000001, 2:6) # no result
#' getIndex(5.00000001, 2:6) # difference is small enough to be found by all.equal
#' @export
getIndex <- function(var, vect){
  index <- which(vect == var) 
  if(length(index) == 0){
    #print("all.equal was used")
    comp <- mapply(function(x) {isTRUE(all.equal(x, var, check.names = FALSE))}, vect) 
    index <- which(comp == TRUE)
  }
  return(index)
}      # get index of variable „var“ in vector „vect“

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
  
  n <- length(obj@init) - 1 # number of continuous variables
  nj <- length(obj@ratesprays) # number of jumptypes
  k <- length(obj@discStates[[1]]) # number of different discrete states
  z <- obj@init
  
  ratematrix <- rep(list(rep(list(NULL), k)), k)
  for(jtype in 1:nj){ 
    for(i in 1:k){
      oldDiscVar <- obj@discStates[[1]][[i]]
      z[length(z)] <- oldDiscVar
      for(j in 1:k){
        newDiscVar <- obj@discStates[[1]][[j]]
        if(obj@jumpfunc(1, z, obj@parms, jtype)[length(z)] == newDiscVar){
          ratematrix[[i]][[j]] <- obj@ratesprays[[jtype]][[i]]
        }
      } 
    }
  }
  return(ratematrix)
}

#' @note only works for one discrete variable
#' @importFrom spray is.zero subs lone
blowupSpray <- function(obj, spray){
  # This function blows the last variable of the spray object 
  # (which stands for the discrete variable θ) up to several different 
  # indicator variables θ₁,…,θₖ, where k = # different states.
  
  # example: x² + 2*θ + θy becomes x² + 0*θ₁ + 2*1*θ₂ + 0*θ₁*y + 1*θ₂*y, if discDomain = c(0,1)
  #          blowupPoly(polyModel9, linear(c(1,0,0),2) + linear(c(0,0,2),1) + product(c(0,1,1)))
  
  if(is.zero(spray)) return(spray)
  
  n <- length(obj@init) # number of continous variables
  k <- length(obj@discStates[[1]]) # number of different discrete states
  spray <- increase_arity(spray, n - 1 + k)
  
  # part of spray that is independent of θ:
  core <- subs(spray, n, 0, keepArity = TRUE) 
  blowedSpray <- core
  
  # add the specific parts (depend on θ):
  for(i in 1:k){
    discVar <- obj@discStates[[1]][i]
    if(!is.zero(spray-core)){
      specific <- subs(spray-core, n, discVar, keepArity = TRUE)*lone(n-1+i, n-1+k)
      blowedSpray <- blowedSpray + specific
    }
  }
  return(blowedSpray)
}  # change discrete variable θ to indicator variables θ₁,…,θₖ


##### output methods ####

setMethod(f="print",
          signature="polyPdmpModel",
          definition=function(x, ...)
          {
            .local <- function (x, all = FALSE, ...) 
            {
              if (all) {
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
            .local(x, ...)
          }
)
