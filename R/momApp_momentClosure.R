#t1 wird die variable mean überhaupt benötigt?

########### moment closure ##############

#' Calculate moments
#'
#' Compute central moments of a multivariate variable \eqn{X = (X_1, ..., X_n)}.
#'
#' @param mean vector or list with expected values. Entry \code{mean[[i]]} should
#' contain a value for \eqn{E(X_i)}.
#' @param cov matrix or nested list. Entry \code{cov[i][j]} (or \code{cov[[i]][[j]]} 
#' respectivly) should contain a #' value for the covariance \eqn{Cov(X_i, X_j)}.
#' @param var optional. If \code{n = 1}, \code{cov} would only
#' contain one entry - the variance \eqn{Var(X)}. This value can be passed
#' to parameter \code{var} instead of \code{cov}.
#' @param missingOrders numeric matrix. Each row gives the order
#' of a moment that shall be calculated.
#' @param distribution string specifying the (multivariate) distribution
#' of X. The following values are possible: \code{zero} sets
#' all moments to 0, \code{normal} sets all moments as moments of
#' a centralized normal distribution.
#' @return A list where each element is a quoted expression.
#' The i-th element of this list gives a formula for the
#' moment whose order is given in the i-th row of \code{missingOrders}.
#' @importFrom symmoments callmultmoments
#' @importFrom stringr str_extract_all
#' @export
momentClosure <- function(distribution, missingOrders, mean = NA, cov = NA, var = NA){
  
  if(is.vector(missingOrders))
    missingOrders <- matrix(missingOrders, nrow = 1)
  n <- ncol(missingOrders)
  
  # validate cov
  if(!is.na(var))
    cov <- var
  if(is.matrix(cov))
    cov <- apply(cov, 1, as.list)
  stopifnot(n == length(cov))
  for(i in 1:n){
    for(j in 1:n)
      if(cov[[i]][[j]] != cov[[j]][[i]])
        stop("The covariance matrix 'cov' should be symmetric.")
  }
  
  # definitions
  missingMoments <- rep(list(NA), nrow(missingOrders))
  n <- ncol(missingOrders)
  
  # zero: all moments = 0
  if(distribution == "zero"){
    missingMoments <- rep(list(0), nrow(missingOrders))
  }
  
  if(distribution == "normal"){
    for(i in seq_len(nrow(missingOrders))){
      order <- as.numeric(missingOrders[i, ])
      if(sum(order) %% 2 == 1){
        missingMoments[[i]] <- 0
      }
      else{
        moment <- symmoments::callmultmoments(order) # calculation of the moment
        names(moment$coefficients) <- NULL
        momFormula <- list()
        for(j in length(moment$coefficients)){
          factors <- lapply(1:(n*(n+1)/2), function(r){
            power <- moment$representation[j, r]
            index <- grep("[0-9]", names(moment$representation[j, ])[r], value = TRUE)
            index <- as.numeric(stringr::str_extract_all(index, "[[:digit:]]")[[1]])
            if(power != 0) bquote(.(cov[[index[1]]][[index[2]]])^.(power))
          })
          factors[sapply(factors, is.null)] <- NULL
          factors <- Reduce(function(a,b) bquote(.(a)*.(b)), factors)
          momFormula <- append(momFormula,
                               bquote(.(moment$coefficients[j])*.(factors)))
        }
        momFormula <- Reduce(function(a,b) bquote(.(a)+.(b)), momFormula)
        missingMoments[[i]] <- momFormula
      }
    }
  }
  
  # normal for one dimensional variable
  if(distribution == "normal1d" & n == 1){
    for(i in seq_len(nrow(missingOrders))){
      order <- as.numeric(missingOrders[i, 1])
      # sigmaRowIndex <- prodlim::row.match(2, knownOrders) #t2 check if 2nd moment is  contained in lhs
      # sigma <- knownMoments[[sigmaRowIndex]]
      variance <- cov[[i]][[1]]
      missingMoments[[i]] <- switch(order %% 2,
                                    bquote(.(variance)^.(order)*.(dfactorial(order-1))),
                                    0)
    }
  }
  
  if(anyNA(missingMoments)){
    stop("Closure method '", closure, "' is not implemented.")
  }
  
  return(missingMoments)
}
