library(spray) # multipol with sparse arrays
options(polyform = TRUE)

#### todo ####
# keepArity in drop umbenennen

##### new version of 'spray' ####
# I changed the handling of the NULL polynomial, 
# because scalar*(the NULL polynomial) didn't work

#' @importFrom spray spraymaker
my.spray <- function (M, x, addrepeats = FALSE){
  if (is.null(M)) { #new
    return(spraymaker(list(NULL, NULL)))
  }
  if (is.vector(M)) {
    M <- rbind(M)
  }
  if (inherits(M, "partition")) {
    M <- t(M)
  }
  M <- as.matrix(M)
  if (missing(x)) {
    x <- rep(1, nrow(M))
  }
  if (length(x) == 1) {
    x <- rep(x, nrow(M))
  }
  return(spraymaker(list(M, x), addrepeats = addrepeats))
}

# assign to "spray"
unlockBinding("spray", as.environment("package:spray"))
assignInNamespace("spray", my.spray, ns="spray", envir=as.environment("package:spray"))
assign("spray", my.spray, as.environment("package:spray"))
lockBinding("spray", as.environment("package:spray"))

##### new version of 'subs' ####
# as method 'subs' didn't work as expected, i had to redefine it

# new version: 
#' @importFrom spray spray value process_dimensions
my.subs <- function (S, dims, x, keepArity = FALSE){
  dims <- process_dimensions(S, dims)
  
  if(keepArity) {
    index <- index(S)
    index[, dims] <- 0
  }
  else index <- index(S)[, -dims, drop = FALSE]
  
  val <- Reduce("*", as.data.frame(x^drop(index(S)[, dims, drop = FALSE])))
  return(spray(index, val * value(S), addrepeats = TRUE))
}

# assign to "subs"
unlockBinding("subs", as.environment("package:spray"))
assignInNamespace("subs", my.subs, ns="spray", envir=as.environment("package:spray"))
assign("subs", my.subs, as.environment("package:spray"))
lockBinding("subs", as.environment("package:spray"))

### old version:
# subs <- function (S, dims, x) 
# {
#   dims <- process_dimensions(S, dims)
#   return(spray(index(S)[, -dims, drop = FALSE], 
#                drop(index(S)[, dims, drop = FALSE] %*% x), 
#                addrepeats = TRUE))
# }

#### increase_arity ####
# increase arity of a spray object

#' @importFrom spray lone spray value
increase_arity <- function(S, arity){
  if(is.null(arity(S))) return(0*lone(1, 1))
  if(arity(S) > arity) stop("The value for 'arity' is to low.")
  index <- cbind(index(S), matrix(0, nrow = length(index(S)[,1]), ncol = arity - arity(S)))
  spray(index, value(S))
}

#### load spray-package #####


