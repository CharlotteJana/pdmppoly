#===============================================================================
#v1 Zur Dokumentation von dtrunc: Zitat aus http://r-pkgs.had.co.nz/check.html
# "If the licenses are compatible you can copy and paste the exported function 
#  into your own package. If you do this, remember to update Authors@R."

#------ dtrunc & mtrunc ----------

#' Probability density function of truncated random variables
#' 
#' This function computes values for the probability density function of a
#' truncated random variable. It was originally implemented in package
#' \pkg{truncdist} and slightly modified to return zeros in case that the
#' trunction interval [a, b] is not inside the support of the density function.
#' @inherit truncdist::dtrunc
#' @param ... other arguments are passed to the corresponding quantile function
#' @examples 
#' x <- seq(0, 3, 0.1)
#' dtrunc(x, spec = "norm", a = 1, b = 2)
#' curve(dtrunc(x, spec = "norm", a = -Inf, b = 1), -10, 2)
#' 
#' \dontrun{
#' # different results for intervals outside the support of the density function:#' 
#' truncdist::dtrunc(x, spec = "norm", a = 20, b = 30) # gives error
#' pdmppoly::dtrunc(x, spec = "norm", a = 20, b = 30) # gives only a warning}
#' @importFrom utils str
#' @export
dtrunc <- function (x, spec, a = -Inf, b = Inf, ...) {
  if (a >= b) 
    stop("argument a is greater than or equal to b")
  tt <- rep(0, length(x))
  g <- get(paste("d", spec, sep = ""), mode = "function")
  G <- get(paste("p", spec, sep = ""), mode = "function")
  if(G(b, ...)-G(a, ...) != 0){ # additional to truncdist::dtrunc
    tt[x >= a & x <= b] <- g(x[x >= a & x <= b], ...)/(G(b, ...) - G(a, ...))
  }
  else  
   warning("Truncation interval is outside the support of the density function")
  
  return(tt)
}

#' Moments of truncated random variables
#' 
#' This function computes the raw moments of a truncated random variable.
#' @inheritParams truncdist::dtrunc
#' @param order numeric vector giving the order of the moments
#' @param ... other arguments are passed to the corresponding moment function
#' @examples
#' mtrunc(1:6, spec = "norm")
#' mtrunc(1:6, spec = "norm", b = 0)
#' @seealso \code{\link{dtrunc}} for the probability distribution of a truncated variable.
#' @importFrom actuar mnorm mlnorm mexp munif
#' @importFrom stats integrate
#' @export
mtrunc <- function (order, spec, a = -Inf, b = Inf , ...){
  argList <- list(...)
  
  if(spec == "norm" & a == -Inf & b == Inf) 
    return(actuar::mnorm(order, ...))
  else if(spec == "lnorm" & a <= 0 & b == Inf)   
    return(actuar::mlnorm(order, ...))
  else if(spec == "exp" & a <= 0 & b == Inf)     
    return(actuar::mexp(order, ...))
  else if(spec == "unif") {
    if(a <= argList$min & b >= argList$max) 
      return(actuar::munif(order, ...))
    if(a > argList$max | b < argList$min) 
      stop("this is not a distribution")
  }
  else{
  #message("Direct integration is used.")
  sapply( order, 
    function(i) integrate( 
      Vectorize(function(x) x^i * dtrunc(x, spec = spec, a = a, b = b, ...)), 
      lower=a, upper=b)$value)
  }
}

#------ dmix & mmix ----------

#' Mixtures of truncated distributions
#' 
#' Compute the distribution (\code{dmix}) and raw moments (\code{mmix}) of a
#' mixture of different, eventually truncated distributions.
#' 
#' @param lower numeric giving the lower bound of the support. Defaults to -Inf.
#' @param upper numeric giving the upper bound of the support. Defaults to Inf.
#' @param distrib a list. Every element is itself a list representing a
#'   distribution. This list should have one named element \code{spec} which
#'   describes the distribution, i. e. "exp" for the exponential distribution or
#'   "norm" for the normal distribution. The other elements are optional
#'   additional parameters for the specfied distribution.
#' @param weights numeric vector with the same length as \code{distrib}.
#' Provides weights for every distribution given in \code{distrib}.
#' @examples 
#' distributions <- list(list(spec="exp", rate = 2),
#'                      list(spec="norm", mean = 0, sd = 0.5),
#'                      list(spec="unif", min = 2, max = 3))
#' d <- dmix(-1, 3, weights = c(.1,.3,.1), distrib = distributions)
#' curve(d(x), -2, 5)
#' @name dmix
#' @aliases mmix
#' @seealso \code{\link{random.distribution}} to plot a randomly created mixture
#' of truncated distributions.
#' @export
dmix <- function(lower = -Inf, upper = Inf, distrib, weights){
  
  n <- length(distrib)
  if(missing(weights)) weights = rep(1/n, n)
  
  weights <- weights/sum(weights)
  function(x){
    h <- sapply(1:n, function(i) 
      weights[i]*do.call("dtrunc", c(x = list(x), distrib[[i]], 
                                     a = lower, b = upper))
    )
    rowSums(h)
  }
}

#' @param order order of the moment
#' @rdname dmix
#' @export
mmix <- function(order, lower = -Inf, upper = Inf, distrib, weights){
  
  n <- length(distrib)
  if(missing(order)) order = 1:4
  if(missing(weights)) weights = rep(1/n, n)
  
  weights <- weights/sum(weights)
  h <- sapply(1:n, function(i) 
    weights[i]*do.call("mtrunc", c(list(order = order), 
                                   a = lower, 
                                   b = upper, 
                                   distrib[[i]]))
  )
  if(length(order) > 1)
    return(rowSums(h))
  else
    return(sum(h))
}

#------ random.distribution ----------

#' Create and plot a random distribution
#' 
#' This function computes a random truncated distribution which is a mixture of
#' one up to 4 different distributions that can be uniform, exponential,
#' normal or lognormal distributed. It can be used to test the evectiveness
#' of function \code{\link{is.unimodal}}.
#' @param lower lower bound of the support of the distribution. Defaults to 0.
#' @param upper upper bound of the support of the distribution. Defaults to 10.
#' @param plot boolean variable. Should the distribution be plotted?
#' @return A list with the distributions that are components of the mixture
#'   distribution. It has the same structure as parameter \code{distrib} in
#'   function \code{\link{dmix}}.
#' @examples
#' par(mfrow = c(1, 2))
#' d <- random.distribution(plot = TRUE)
#' curve(dmix(0, 10, distrib = d)(x), -1, 10, type = "l", col = "red")
#' @importFrom stats runif
#' @export
random.distribution <- function(lower = 0, upper = 10, plot = TRUE){
  x <- NULL # to avoid warning 'no visible binding...' in R CMD check
  n <- sample(1:4, 1)                         # random number of components
  specs <- c("norm", "unif", "exp", "lnorm")  # possible shape of components
  weights <- sample(1:3, n, replace=TRUE)     # random weights
  
  distributions <- as.list(NULL) # random values for distribution parameters
  for(i in 1:n){
    unifMin <- sample(lower:((lower+upper)/2),1)
    t <- switch(sample(specs,1),
                "norm"  = list(spec="norm", mean = sample(i:(2*i),1)),
                "unif"  = list(spec="unif", min = unifMin, max = unifMin+1),
                "exp"   = list(spec="exp", rate = runif(1, min = 0.1, max = 2)),
                "lnorm" = list(spec="lnorm", meanlog = sample(i:(2*i),1)))
    distributions[[i]] <- t
  }
  if(plot) graphics::curve(dmix(distrib=distributions, lower = lower, upper = upper)(x), 
                 from = lower-1, to = upper+1, ylab = "distribution")
  return(distributions)
}