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
#' @value A list where each element is a quoted expression.
#' The i-th element of this list gives a formula for the
#' moment whose order is given in the i-th row of \code{missingOrders}.
momentClosure <- function(distribution, missingOrders, mean = NA, cov = NA, var = NA){
  
  if(is.vector(missingOrders))
    missingOrders <- matrix(missingOrders, nrow = 1)
  if(!is.na(var))
    cov <- var
  if(is.matrix(cov))
    cov <- apply(cov, 1, as.list)
  
  # definitions
  missingMoments <- rep(list(NA), nrow(missingOrders))
  n <- ncol(missingOrders)
  
  # zero: all moments = 0
  if(distribution == "zero"){
    missingMoments <- rep(list(0), nrow(missingOrders))
  }
  
  # normal for one dimensional variable
  if(distribution == "normal" & n == 1){
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
  
  # normal: moments of a centralized normal distribution
  if(distribution == "normal"){
    for(i in seq_len(nrow(missingOrders))){
      order <- missingOrders[i, ]
      if(sum(order) %% 2 == 1){
        missingMoments[[i]] <- 0
      }
      else{
       k <- sum(order)
       pairs <- combn(1:k, 2)
       pairCombinations <- combn(1:ncol(pairs), k/2)
       combinations <- NULL
       for(i in seq_len(ncol(pairCombinations))){
         combi <- as.vector(sapply(pairCombinations[, i], function(j) pairs[, j]))
         if(length(unique(combi)) == length(combi))
           combinations <- rbind(combinations, combi)
       }
       I <- rep(1:n, order)
       combinations <- apply(combinations, 2, function(i) I[i])
       sum <- list()
       for(i in seq_len(nrow(combinations))){
         factor <- lapply(1:(k/2), function(j) cov[[combinations[i, 2*j-1]]][[combinations[i,2*j]]])
         factor <- Reduce(function(a,b) bquote(.(a)*.(b)), factor)
         sum <- append(sum, factor)
       }
       sum <- Reduce(function(a,b) bquote(.(a)+.(b)), sum)
       return(sum)
      }
    }
  }
  
  if(anyNA(missingMoments)){
    stop("Closure method '", closure, "' is not implemented.")
  }
  
  return(missingMoments)
}
