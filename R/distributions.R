#======== todo =================================================================
#t1 documentation für alle funktionen
#t1 dmix: parameter distrib beschreiben
#t1 wie dtrunc zitieren?
#t2 älteren Code rausnehmen -> wohin?

#library(truncdist) # enthält die Funktion dtrunc

###### Momente für beliebige gestutzte Verteilungen ######

#' @importFrom utils str
#' @export
dtrunc <- function (x, spec, a = -Inf, b = Inf, ...) {
  # dtrunc stammt aus truncdist und wurde nur geringfügig von mir geändert
  if (a >= b) 
    stop("argument a is greater than or equal to b")
  tt <- rep(0, length(x))
  g <- get(paste("d", spec, sep = ""), mode = "function")
  G <- get(paste("p", spec, sep = ""), mode = "function")
  if(G(b, ...)-G(a, ...) != 0){ # zusätzlich zu dtrunc
    tt[x >= a & x <= b] <- g(x[x >= a & x <= b], ...)/(G(b, ...) - G(a, ...))
  }
  else {print(paste(spec, "a =", a, "b =", b)); str(list(...)); 
                    stop("this is not a distribution")}
  return(tt)
}

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

#------ dmix ----------

#' Mixtures of truncated distributions
#' 
#' @param a numeric giving the lower bound of the support. Defaults to -Inf.
#' @param b numeric giving the upper bound of the support. Defaults to Inf.
#' @param distrib
#' @param weights numeric vector with the same length as \code{distrib}.
#' Provides weights for every distribution given in \code{distrib}.
#' @examples 
#' distributions <- list(list(spec="exp", rate = 2),
#'                      list(spec="norm", mean=0, sd = 0.5),
#'                      list(spec="unif", min=2, max=3))
#' curve(dmix(3, a=-1, weights = c(.1,.3,.1), distrib=distributions)(x),-2,5)
#' @name dmix
#' @aliases mmix
#' @export
dmix <- function(a = -Inf, b = Inf, distrib, weights){
  
  n <- length(distrib)
  if(missing(weights)) weights = rep(1/n, n)
  
  weights <- weights/sum(weights)
  function(y){
    h <- sapply(1:n, function(i) 
      weights[i]*do.call("dtrunc", c(x=list(y), a=a, b=b, distrib[[i]]))
    )
    rowSums(h)
  }
}

#' @rdname dmix
#' @export
mmix <- function(order, a = -Inf, b = Inf, weights, distrib){
  
  n <- length(distrib)
  if(missing(order)) order = 1:4
  if(missing(weights)) weights = rep(1/n, n)
  
  weights <- weights/sum(weights)
  h <- sapply(1:n, function(i) 
    weights[i]*do.call("mtrunc", c(list(order=order), a=a, b=b, distrib[[i]]))
  )
  if(length(order) > 1) 
    h <- rowSums(h)
  else 
    h <- sum(h)
  return(h)
}

#' @importFrom stats runif
#' @export
random.distribution <- function(A = 0, B = 10, curve = TRUE){
  x <- NULL # to avoid warning 'no visible binding...' in R CMD check
  n <- sample(1:4, 1)                         # random number of components
  specs <- c("norm", "unif", "exp", "lnorm")  # possible shape of components
  weights <- sample(1:3, n, replace=TRUE)     # random weights
  
  distributions <- as.list(NULL) # random values for distribution parameters
  for(i in 1:n){
    unifMin <- sample(A:((A+B)/2),1)
    t <- switch(sample(specs,1),
                "norm"  = list(spec="norm", mean = sample(i:(2*i),1)),
                "unif"  = list(spec="unif", min = unifMin, max = unifMin+1),
                "exp"   = list(spec="exp", rate = runif(1, min = 0.1, max = 2)),
                "lnorm" = list(spec="lnorm", meanlog = sample(i:(2*i),1)))
    distributions[[i]] <- t
  }
  if(curve) curve(dmix(distrib=distributions, a = A, b = B)(x), A-1, B+1)
  return(distributions)
}

##### Kurzformen für Mixtures gleichen Typs #######

# #Normalverteilungen
# dmixnorm <- function(n, dist=5, means, sds, ...){
#   if(missing(means)) means = dist*1:n
#   if(missing(sds)) sds = rep(1, n)
#   distributions <- lapply(1:n, function(i) list(spec = "norm", mean = means[i], sd = sds[i]))
#   function(x) dmix(distrib = distributions, ...)(x)
# }
# mmixnorm <- function(n, order, dist=5, means, sds, ...){
#   if(missing(means)) means = dist*1:n
#   if(missing(sds)) sds = rep(1, n)
#   distributions <- lapply(1:n, function(i) list(spec = "norm", mean = means[i], sd = sds[i]))
#   mmix(order, distrib = distributions, ...)
# }
# 
# #Gleichverteilungen
# dmixunif <- function(n, dist = 1, length = 1, lowers, uppers, ...){
#   if(missing(lowers)) lowers = (dist+length)*1:n
#   if(missing(uppers)) uppers = (dist+length)*1:n + length
#   distributions <- lapply(1:n, function(i) list(spec = "unif", min = lowers[i], max = uppers[i]))
#   function(x) dmix(distrib = distributions, ...)(x)
# }
# mmixunif <- function(n, order, dist = 1, length = 1, lowers, uppers, ...){
#   if(missing(lowers)) lowers = (dist+length)*1:n
#   if(missing(uppers)) uppers = (dist+length)*1:n + length
#   distributions <- lapply(1:n, function(i) list(spec = "unif", min = lowers[i], max = uppers[i]))
#   mmix(order, distrib = distributions, ...)
# }
# 
# #### LogNormalverteilungen ####  
# dmixlnorm <- function(n, means, sds, ...){
#   if(missing(means)) means = 1:n
#   if(missing(sds)) sds = rep(c(1/3, 3), n)
#   distributions <- lapply(1:n, function(i) list(spec = "lnorm", meanlog = means[i], sdlog = sds[i]))
#   function(x) dmix(distrib = distributions, ...)(x)
# }
# mmixlnorm <- function(n, order, means, sds, ...){
#   if(missing(means)) means = 1:n
#   if(missing(sds)) sds = rep(c(1/3, 3), n)
#   distributions <- lapply(1:n, function(i) list(spec = "lnorm", meanlog = means[i], sdlog = sds[i]))
#   mmix(order, distrib = distributions, ...)
# }


##### gestutzte Normalverteilungennach Orjebin #####

#Momente der gestutzten Normalverteilung nach Orjebin_2014

#' @importFrom stats dnorm qnorm pnorm
#' @export
mtnorm <- function(order, mean = 0, sd = 1, lower = - Inf, upper = Inf){
  max <- max(order)
  moments <- 0:(max+1) # startet mit m₋₁ = 0, and m₀ = 1
  k <- 3
  d <- function(x) dnorm(x, mean = mean, sd = sd)
  p <- function(x) pnorm(x, mean = mean, sd = sd)
  uppern <- (upper-mean)/sd
  lowern <- (lower-mean)/sd
  while(k <= max+2){
   moments[k] <- (k-3)*sd^2*moments[k-2]+mean*moments[k-1] - 
                sd*(upper^(k-3)*d(upper)-lower^(k-3)*d(lower))/
                (p(upper)-p(lower))
   k <- k+1
  }
  print(d(-Inf)*(-Inf))
  moments[order+2]
}

#seltsames Verhalten für mean, bsp:  
#mtnorm(1:15, mean = 20, lower = -600, upper = 600)-mnorm(1:15, mean = 20)
#gibt große Differenz, die auch bei weiteren Grenzen nicht verschwindet