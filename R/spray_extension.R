library(spray) # multipol with sparse arrays
options(polyform = TRUE)

#### todo ####
# keepArity in drop umbenennen

##### new version of 'spray' ####
# I changed the handling of the NULL polynomial, 
# because scalar*(the NULL polynomial) didn't work

#' #' @importFrom spray spraymaker
#' my.spray <- function (M, x, addrepeats = FALSE){
#'   if (is.null(M)) { #new
#'     return(spraymaker(list(NULL, NULL)))
#'   }
#'   if (is.vector(M)) {
#'     M <- rbind(M)
#'   }
#'   if (inherits(M, "partition")) {
#'     M <- t(M)
#'   }
#'   M <- as.matrix(M)
#'   if (missing(x)) {
#'     x <- rep(1, nrow(M))
#'   }
#'   if (length(x) == 1) {
#'     x <- rep(x, nrow(M))
#'   }
#'   return(spraymaker(list(M, x), addrepeats = addrepeats))
#' }
#' 
#' # assign to "spray"
#' unlockBinding("spray", as.environment("package:spray"))
#' assignInNamespace("spray", my.spray, ns="spray", envir=as.environment("package:spray"))
#' assign("spray", my.spray, as.environment("package:spray"))
#' lockBinding("spray", as.environment("package:spray"))

##### new version of 'subs' ####
# as method 'subs' didn't work as expected, i had to redefine it

# new version: 
#' #' @importFrom spray spray value process_dimensions
#' my.subs <- function (S, dims, x, keepArity = FALSE){
#'   dims <- process_dimensions(S, dims)
#'   
#'   if(keepArity) {
#'     index <- index(S)
#'     index[, dims] <- 0
#'   }
#'   else index <- index(S)[, -dims, drop = FALSE]
#'   
#'   val <- Reduce("*", as.data.frame(x^drop(index(S)[, dims, drop = FALSE])))
#'   return(spray(index, val * value(S), addrepeats = TRUE))
#' }
#' 
#' # assign to "subs"
#' unlockBinding("subs", as.environment("package:spray"))
#' assignInNamespace("subs", my.subs, ns="spray", envir=as.environment("package:spray"))
#' assign("subs", my.subs, as.environment("package:spray"))
#' lockBinding("subs", as.environment("package:spray"))

### old version:
# subs <- function (S, dims, x) 
# {
#   dims <- process_dimensions(S, dims)
#   return(spray(index(S)[, -dims, drop = FALSE], 
#                drop(index(S)[, dims, drop = FALSE] %*% x), 
#                addrepeats = TRUE))
# }

#' Increase arity of a spray object
#'
#' Increase the arity of a spray object by adding new columns (with zeros) to
#' the index matrix of the spray object.
#'
#' The arity of a spray object defines the dimension of the polynomial. For
#' example \code{lone(1, 2)} and \code{lone(1, 5)} look the same when printed to
#' the console, but differ in their arity (2 or 5). This can be seen by using
#' method \code{\link[spray]{print_spray_matrixform}}. Calculation with spray
#' objects (sum, product etc.) is only possible if the spray objects havte the
#' same arity. This is where function \code{increase_arity} comes into play. It
#' can increase the arity of a polynomial/spray object by adding new variables.
#' This means essentially inserting new columns with zeros in the index matrix
#' of the spray object. Parameter \code{additionalCols} defines the columnindex
#' of the these additinoal columns.
#' @param S spray object that shall have a new arity
#' @param additionalCols column indexes where to add new variables
#' @examples
#' lone(1, 5) + lone(1, 2) # throws error
#' lone(1, 5) + increase_arity(lone(1, 2), 3:5)
#' 
#' print_spray_matrixform(homog(2, 3))
#' s <- increase_arity(homog(2, 3), 2) # insert variable in the middle
#' print_spray_matrixform(s)
#' @importFrom spray lone spray value
increase_arity <- function(S, additionalCols){
  
  n <- length(additionalCols)
  if(!all(additionalCols %in% 1:(arity(S) + n))) 
    stop("The values of additionalCols do not match with the arity of S.")
  
  arity <- n + ncol(index(S))
  if(is.null(arity(S))) return(0*lone(1, 1)) # oder length(additionalCols)

  
  matrix <- cbind(index(S), matrix(0, nrow = length(index(S)[,1]), ncol = n))
  print(matrix)
  id <- c(1:arity(S), additionalCols - 0.5)
  print(id)
  matrix <- matrix[, order(id)]
  spray(matrix, value(S))
}

# increase_arity(subs(S, dims, x), dims)