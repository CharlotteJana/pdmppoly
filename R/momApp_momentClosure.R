########### moment closure ##############

#' @param knownOrders numeric matrix. Each row gives the order
#' of a moment that is already known (i. e. moments of order
#' one or two which are needed for calculation).
#' @param knownMoments list with moments. The i-th entry contains
#' the moment of order \code{knownOrders[i, ]}.
#' @param missingOrders numeric matrix. Each row gives the order
#' of a moment that shall be calculated.
#' @param closure string specifying the closure method that is 
#' to be used. The following values are possible: \code{zero} sets
#' all moments to 0, \code{normal} sets all moments as moments of
#' a centralized normal distribution.
#' @value A list where each element is a quoted expression.
#' The i-th element of this list gives a formula for the
#' moment whose order is given in the i-th row of \code{missingOrders}.
momentClosure <- function(closure, missingOrders, knownOrders, knownMoments){
  missingMoments <- rep(list(NA), nrow(missingOrders))
  
  # zero: all moments = 0
  if(closure == "zero"){
    missingMoments <- rep(list(0), length(missingMoments))
  }
  
  # normal: moments of a centralized normal distribution
  n <- ncol(missingOrders)
  if(closure == "normal" & n == 1){ # if we have only one continous variable
    for(i in seq_len(nrow(missingOrders))){
      order <- as.numeric(missingOrders[i, 1])
      sigmaRowIndex <- prodlim::row.match(2, knownOrders) #t2 check if 2nd moment is  contained in lhs
      sigma <- knownMoments[[sigmaRowIndex]]
      missingMoments[[i]] <- switch(order %% 2,
                                    bquote(.(sigma)^.(order)*.(dfactorial(order-1))),
                                    0)
    }
  }
  
  if(anyNA(missingMoments)){
    stop("Closure method '", closure, "' is not implemented.")
  }
  
  return(missingMoments)
}
