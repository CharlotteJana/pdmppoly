#' Double factorial
#' 
#' This function returns the double factorial n!! of a natural number n.
#' It is defined as \eqn{n!! = n*(n-2)*(n-4)*...*1} if n is odd and
#' \eqn{n!! = n*(n-2)*(n-4)*...*2} if n is even.
#' 
#' @param n a non negative, natural number
#' @export
dfactorial <- function(n){
  if(n %% 1 != 0 | n < 0)
    stop("dfactorial is only defined for non negative, natural numbers")
  if(n == 0)
    return(1)
  if(n %% 2 == 1){
    odds <- seq_len(n)[as.logical(seq_len(n) %% 2)]
    return(Reduce("*", odds))
  }
  if(n %% 2 == 0){
    evens <- seq_len(n)[!as.logical(seq_len(n) %% 2)]
    return(Reduce("*", evens))
  }
}

#' Multinomial coefficient
#' 
#' This function returns the multinomial coefficient of a natural number m
#' and a vector of natrual numbers k.
#' It is defined as \eqn{m!\(k[1]!*...*k[n]!)}, 
#' where \code{n = length(k)} and \eqn{m!} stands for \code{factorial(m)}.
#' 
#' @param m a non negative, natural number
#' @param k a vector of non negative, natural numbers
#' @export
multinomial <- function(m, k){
  a <- sapply(k, factorial)
  return(factorial(m)/prod(a))
}
