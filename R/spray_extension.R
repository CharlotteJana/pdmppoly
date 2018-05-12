# library(spray) # multipol with sparse arrays
options(polyform = TRUE)


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
#' @param position column indexes where to add new variables
#' @examples
#' lone(1, 5) + lone(1, 2) # throws error
#' lone(1, 5) + increase_arity(lone(1, 2), 3:5)
#' 
#' print_spray_matrixform(homog(2, 3))
#' s <- increase_arity(homog(2, 3), 2) # insert variable in the middle
#' print_spray_matrixform(s)
#' 
#' p <- subs(product(1:4), dims = 2:3, x = c(-1, 2)) # substitute values
#' increase_arity(p, 2:3) # keep arity as before the substitution
#' @importFrom spray lone spray value
increase_arity <- function(S, position){

  n <- length(position)
  newArity <- n + arity(S)
  
  if(!all(position %in% 1:newArity)) 
    stop("The values of position do not match with the arity of S.")
  
  if(is.null(arity(S))) return(0*lone(1, max(1, newArity)))
  
  Scols <- index(S)
  matrix <- NULL
  for(i in 1:newArity){
    if(i %in% position) # add column with zeros
      matrix <- cbind(matrix, rep(0, nrow(index(S))))
    else{ # add column from polynomial S
      matrix <- cbind(matrix, Scols[, 1])
      if(!is.null(Scols)) 
        Scols <- as.matrix(Scols[, -1], nrow = nrow(index(S)))
    }
  }
  
  
  # br <- which(diff(position) > 1)
  # print(br)
  # p <- cut(position, 
  #          breaks = c(0, position[unique(c(br, length(position)))]), 
  #          right = TRUE)
  # posList <- split(position, p)
  # print(posList)
  # 
  # id <- 1:arity(S)
  # for(i in seq_along(posList)){
  #   k <- length(posList[[i]])
  #   id <- c(id, posList[[i]][1] - k:1/(k+1))
  # }
  # print(id)
  # matrix <- cbind(index(S), matrix(0, nrow = length(index(S)[,1]), ncol = n))
  # matrix <- matrix[, order(id)]
  
  # names(position) <- rep("N", n)
  # idx <- 1:arity(S)
  # names(idx) <- rep("S", arity(S))
  # 
  # # m <- t(index(S))
  # r <- nrow(index(S))
  # for(i in seq_along(position)){ # add columns with zeros at position
  #   idx <- append(idx, position[i], after = position[i]-1)
  #   #m <- matrix(unlist(append(m, rep(0, r), after = position[i]-1)), ncol = r)
  #   
  # }
  spray(matrix, value(S))
}