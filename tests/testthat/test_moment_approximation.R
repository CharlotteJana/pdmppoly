#======== todo =================================================================
#t1 ersten test als demo und dann vereinfachen
#t1 erster test: l > 4 setzen und richtige ergebnisse erhalten


context("moment approximation")

test_that("moment calculation works for model 1", {
  
  ### definitions
  data(genePoly1)
  times(genePoly1) <- c(from = 0, to = 2000, by = 1)
  states <- discStates(genePoly1)[[1]]
  l <- 4 # geht für höhere l werte nicht mehr
  k <- length(states)
  n <- length(genePoly1@init) - 1
  
  ### moment approximation
  momApp <- momentApprox(genePoly1, l)
  last <- nrow(momApp$contRes)
  
  ### create matrix with all moment combinations we are interested in
  s <- NULL
  for(i in 1:n){ 
    m <-  matrix(data = 0, nrow = l, ncol = n+3)
    m[,i] <- 1:l
    s <- rbind(s, m)
  }
  s <- rbind(s, c(rep(0, n), 1))
  colnames(s) <- c(names(genePoly1@init), "calc", "approx")
  s <- as.data.frame(s)
  
  ### continous variables
  for(i in 1:(n*l)){
    m <- s[i,1]
    # theoretical values:
    s[i, n+2] <- with(as.list(genePoly1@parms),
      (α/β)^m*Reduce("*", sapply(0:(m-1), function(j) (κ01+j*β)/(κ01+κ10+j*β))))
    # calculated values:
    s[i, n+3] <- momApp$contRes[last, 
                                prodlim::row.match(as.vector(s[i,1:n]), 
                                                   as.matrix(momApp$contInd))]
  }
  
  ### discrete variables
  # theoretical values:
  s[n*l+1, n+2] <- with(as.list(genePoly1@parms), κ01/(κ01+κ10))
  # calculated values:
  s[n*l+1, n+3] <- states %*% momApp$discRes[last, 2:ncol(momApp$discRes)]
 
  #print(s, digits = 4)
  for(i in 1:nrow(s)){
    expect_equal(s[i, 3], s[i, 4], tolerance = 1e-04)
  }
})

##### compareMomentApproximation (only Model 1) #####

compareMomApp <- function(momApp){
  # momApp = result of momentApprox(polyModel1, l = ...) (has class "momApp")
  # the approximation results for polyModel1 are calculated with a formula
  # (only possible for model1) and compared to the results stored in "momApp"
  
  ### definitions
  model <- polyModel1
  l <- momApp$degree
  k <- ncol(momApp$discRes) - 1
  n <- length(model@init) - 1
  
  ### create matrix with all moment combinations we are interested in
  s <- NULL
  for(i in 1:n){ 
    m <-  matrix(data = 0, nrow = l, ncol = n+3)
    m[,i] <- 1:l
    s <- rbind(s, m)
  }
  s <- rbind(s, c(rep(0, n), 1))
  colnames(s) <- c(names(model@init), "calc", "approx")
  s <- as.data.frame(s)
  
  ### fill calculated moments
  for(i in 1:(n*l)){ # continuous variables
    m <- s[i,1]
    s[i, n+2] <- with(as.list(model@parms),
                      (α/β)^m*Reduce("*", sapply(0:(m-1), function(j) (κ01+j*β)/(κ01+κ10+j*β) ))
    )}
  s[n*l+1, n+2] <- with(as.list(model@parms), κ01/(κ01+κ10)) # discrete variable
  
  ### fill approximated moments
  for(i in 1:(n*l)){ # continuous variables
    s[i, n+3] <- momApp$contRes[nrow(momApp$contRes), row.match(as.vector(s[i,1:n]), as.matrix(momApp$contInd))]
  }
  s[n*l+1, n+3] <- model@discDomain %*% momApp$discRes[nrow(momApp$discRes), 2:ncol(momApp$discRes)] # discrete variable
  
  print(s)
}
