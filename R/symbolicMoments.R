#t1 tests
#t2 lognormal: überprüfen, ob mean wohldefiniert
#v1 was passiert bei falschen werten für mean, cov bei gamma?
#t1 References: LakatosCo2015
#v1 Ist es konsistenter, wenn quote(0) anstatt 0 zurückgegeben wird?
#t1 multinomial testen

#' Symbolic calculation of moments
#'
#' Compute the moments of a multivariate variable \eqn{X = (X_1, ..., X_n)}
#' and return the formula as quoted expression.
#'
#' @param mean vector or list with expected values. Entry \code{mean[[i]]} should
#' contain a value for \eqn{E(X_i)}.
#' @param cov matrix or nested list. Entry \code{cov[i][j]} (or \code{cov[[i]][[j]]} 
#' respectivly) should contain a value for the covariance \eqn{Cov(X_i, X_j)}.
#' @param var optional. If \code{n = 1}, \code{cov} would only
#' contain one entry - the variance \eqn{Var(X)}. This value can be passed
#' to parameter \code{var} instead of \code{cov}.
#' @param missingOrders numeric vector or matrix. Each row gives the order
#' of a moment that shall be calculated.
#' @param distribution string specifying the (multivariate) distribution
#' of X. The following values are possible: 
#' \itemize{
#' \item \code{"zero"} sets all moments to 0,
#' \item \code{"normal"} calculates the moments of a centralized multivariate normal distribution,
#' \item \code{"lognormal"} calculates the raw moments of a multivariate lognormal distribution,
#' \item \code{"gamma"} calculates the raw moments of a multivariate gamma distribution.
#' } 
#' @return A list where each element is a quoted expression.
#' The i-th element of this list gives a formula for the
#' moment whose order is given in the i-th row of \code{missingOrders}.
#' @note The calculation of the central moments of a multivariate normal distribution
#' is based on function \code{\link[symmoments]{callmultmoments}} of package \pkg{symmoments}.
#' If \code{distribution = 'gamma'} and the evaluation of the results leads to NaNs, 
#' then the values of cov and mean do not fit to a multivariate gamma distribution.
#' All entries should of cov should be positive and the diagonals should be large
#' with respect to the other entries. More specifically, the inequations 
#' \eqn{mean[i] > \sum_{k \neq i} \frac{mean[k]*cov[i,k]}{cov[i,i]}}
#' should be satisfied for all i in 1:n.
#' @importFrom symmoments callmultmoments
#' @importFrom stringr str_extract_all
#' @export
symbolicMoments <- function(distribution, missingOrders, mean = NA, cov = NA, var = NA){
  
  # definitions
  if(is.vector(missingOrders))
    missingOrders <- matrix(missingOrders, nrow = 1)
  n <- ncol(missingOrders)
  missingMoments <- rep(list(NA), nrow(missingOrders))
  
  # zero: all moments = 0
  if(distribution == "zero"){
    missingMoments <- rep(list(0), nrow(missingOrders))
    return(missingMoments)
  }
  
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
  
  
  # normal: moments are 0 (if sum(order) is odd) or calculated with package symmoments (if sum(order) is even)
  if(distribution == "normal"){
    for(i in seq_len(nrow(missingOrders))){
      order <- as.numeric(missingOrders[i, ])
      if(sum(order) %% 2 == 1){
        missingMoments[[i]] <- 0
      }
      else{
        moment <- symmoments::callmultmoments(order) # calculation of the moment
        names(moment$coefficients) <- NULL
        cnames <- colnames(moment$representation)
        cindexes <- lapply(stringr::str_extract_all(cnames, "[[:digit:]]"), as.numeric)
        momFormula <- list()
        for(j in seq_along(moment$coefficients)){
          factors <- lapply(1:(n*(n+1)/2), function(r){
            power <- moment$representation[j, r]
            if(power != 0) bquote(.(cov[[cindexes[[r]][1]]][[cindexes[[r]][2]]])^.(power))
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
    return(missingMoments)
  }
  
  # lognormal
  if(distribution == "lognormal"){ 
    sigma <- cov
    mu <- list()
    for(i in 1:n){
      sigma[[i]][[i]] <- list(bquote(log(1+.(cov[[i]][[i]])/.(mean[i])^2)))
      sigma[[i]][[i]] <- sigma[[i]][[i]][[1]]
      mu[[i]] <- bquote(log(.(mean[i]))-0.5*.(sigma[[i]][[i]]))
      for(j in seq_len(i-1)){
        sigma[[i]][[j]] <- bquote(
          log(1 + .(cov[[i]][[j]])/exp(.(mu[[i]])+.(mu[[j]])+0.5*(.(sigma[[i]][[i]])*.(sigma[[j]][[j]]))))
          )
        sigma[[j]][[i]] <- sigma[[i]][[j]]
      }
    }
    for(i in seq_len(nrow(missingOrders))){
      order <- as.numeric(missingOrders[i, ])
      summand1 <- lapply(1:n, function(i) bquote(.(order[i])*.(mu[[i]])))
      summand1 <- Reduce(function(a,b) bquote(.(a)+.(b)), summand1)
      summand2 <- lapply(1:n, function(i) lapply(i:n, function(j){ 
                                                        if(i == j) 
                                                          bquote(.(order[i])^2*.(sigma[[i]][[i]]))
                                                        else 
                                                          bquote(2*.(order[i])*.(order[j])*.(sigma[[i]][[j]]))
                                                 })
                         )
      summand2 <- Reduce(c, summand2) # flatten the nested list
      summand2 <- Reduce(function(a,b) bquote(.(a)+.(b)), summand2)
      missingMoments[[i]] <- bquote(exp(.(summand1)+0.5*.(summand2)))
    }
  }
  
  # gamma
  if(distribution == "gamma"){
    
    # define parameters beta and A
    beta <- lapply(1:n, function(i) bquote(.(cov[[i]][[i]])/.(mean[[i]])))
    A_indexes <- as.data.frame(t(cbind(combn(1:n, 2), matrix(rep(1:n, each = 2), nrow = 2))))
    A <- list()
    for(r in 1:nrow(A_indexes)){
      i <- A_indexes[r, 1]
      j <- A_indexes[r, 2]
      if(i != j){
        A <- append(A, bquote(.(cov[[i]][[j]])/(.(beta[[i]])*.(beta[[j]]))))
      }
      else{
        sum_indexes <- which(rowSums(apply(A_indexes, 2, match, i, nomatch = 0)) == 1)
        sum <- Reduce(function(a,b) bquote(.(a)+.(b)), A[sum_indexes])
        A <- append(A, bquote(.(mean[[i]])^2/.(cov[[i]][[i]]) - .(sum)))
      }
    }
    
    # calculate moments
    for(i in seq_len(nrow(missingOrders))){
      order <- as.numeric(missingOrders[i, ])
      polynomials <- list()
      for(j in 1:n){
        poly_indexes <- rowSums(apply(A_indexes, 2, match, j, nomatch = 0)) >= 1
        polynomials <- append(polynomials, list(linear(as.numeric(poly_indexes), 1)^order[j]))
      }
      polynomials <- Reduce('*', polynomials)
      summands <- list()
      for(r in seq_len(nrow(polynomials$index))){
        row <- polynomials$index[r, ]
        factors <- list(polynomials$value[r])
        for(k in seq_along(row)){
          if(row[k] == 1)
            factors <- append(factors, bquote(.(A[[k]])))
          if(row[k] > 1)
            factors <- append(factors, bquote(factorial(.(A[[k]])+.(row[k])-1)/factorial(.(A[[k]])-1)))
        }
        summands <- append(summands, Reduce(function(a,b) bquote(.(a)*.(b)), factors))
      }
      sum <-  Reduce(function(a,b) bquote(.(a)+.(b)), summands)
      result <- append(sum, lapply(seq_along(beta), function(j) bquote(.(beta[[j]])^.(order[j])))) # times beta^order
      missingMoments[[i]] <- Reduce(function(a,b) bquote(.(a)*.(b)), result)
    }
  }
  
  # normal for one dimensional variable
  if(distribution == "normal1d" & n == 1){
    for(i in seq_len(nrow(missingOrders))){
      order <- as.numeric(missingOrders[i, 1])
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
