#======== todo =================================================================
#t3 options(polyform = TRUE) setzen?

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
#' 
#' q <- subs(lone(1,1), 1, 0) # NULL polynomial without arity
#' increase_arity(q, c(2, 5))
#' @note If \code{S} is the NULL polynomial without a given arity (this can
#' arise for example by using method \code{\link[spray]{subs}}), the returned
#' spray object will be the NULL polynomial with its arity given as the maximal
#' entry in vector \code{positions}. 
#' @importFrom spray lone spray value
increase_arity <- function(S, position){

  n <- length(position)
  newArity <- ifelse(is.null(arity(S)),
                     max(position),
                     n + arity(S))
  
  if(!all(position %in% 1:newArity)) 
    stop("The values of position do not match with the arity of S.")
  
  if(is.null(arity(S))) return(0*lone(1, max(1, newArity)))
  
  matrix <- NULL
  j <- 1
  for(i in 1:newArity){
    if(i %in% position) # add column with zeros
      matrix <- cbind(matrix, rep(0, nrow(index(S))))
    else{ # add column from index matrix from S
      matrix <- cbind(matrix, index(S)[, j])
      j <- j+1
    }
  }

  spray(matrix, value(S))
}