#======== todo =================================================================

compute.central.moment <- function(j, moments){
  e <- moments[1]
  s <- sapply(1:j, function(k) choose(j,k)*moments[k]*(-e)^(j-k))
  sum(s) + (-e)^j
}

centr.moments <- function(moments){
  sapply(1:length(moments), function(j) compute.central.moment(j, moments))
}

stand.moments <- function(moments){
  sd <- sqrt(moments[2]-moments[1]^2)
  sapply(1:length(moments), function(j) compute.central.moment(j,moments)/sd^j)
}



random.distribution <- function(A=0, B=10, curve=TRUE){
  n <- sample(1:4, 1)                         # zufällige Anzahl an Komponenten ϵ {1,2}
  specs <- c("norm", "unif", "exp", "lnorm")  # mögliche Gestalt der Komponenten
  weights <- sample(1:3, n, replace=TRUE)     # zufällige Gewichte
  
  distributions <- as.list(NULL)              # weitere zufällige Kenngrößen der Komponenten
  for(i in 1:n){
    unifMin <- sample(A:((A+B)/2),1)
    t <- switch(sample(specs,1),
                "norm"  = list(spec="norm", mean = sample(i:(2*i),1)),
                "unif"  = list(spec="unif", min = unifMin, max = unifMin+1),
                "exp"   = list(spec="exp", rate = runif(1, min = 0.1, max = 2)),
                "lnorm" = list(spec="lnorm", meanlog = sample(i:(2*i),1)))
    distributions[[i]] <- t
  }
  if(curve) curve(dmix(distrib=distributions, a=A, b=B)(x), A-1, B+1)
  return(distributions)
}

test.momentboundaries <- function(n = 100){
  A <- 0                                      # untere Grenze = 0
  B <- sample(1:50,1)                         # zufällige obere Grenze ϵ {1, 2, ..., 50}
  test <- function(i){
    m <- mmix(1:6, A, B, distrib = random.distribution(A, B, curve=FALSE))
    if(!compare.momentboundaries(A, B, m)){
      str(distributions)
      print(m)
      print(paste("B =",B))
      print(paste("Integral:", integrate(dmix(distrib=distributions, a=A, b=B), lower = A, upper = B)$value))
      return(FALSE)
    }
    return(TRUE)
  }
  sapply(1:n, test)
}