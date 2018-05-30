#======== todo =================================================================
#t3 delete old code

library(spray)
context("generator")

#------------ helper functions ---------------

funcdiff <- function(f, val1, val2){
  increase_arity(subs(f, arity(f), val1), arity(f)) - 
    increase_arity(subs(f, arity(f), val2), arity(f))
}

printGenerators <- function(polyMod, formalGen, m){
  
  states <- polyMod@discStates[[1]]
  n <- length(polyMod@init) - length(polyMod@discStates)
  k <- length(states)
  
  polyGen <- function(discVar)
    polyGenerator(polyMod)(product(c(m, 0)))(discVar)
  EVGen <- function(discVar) 
    EVGenerator(polyMod, m, discVar)
  
  for(i in 1:k){
    cat(noquote("polyGenerator    \tdiscVar ="), states[i], "\n")
    print(polyGen(states[i]))
    cat(noquote("formal generator \tdiscVar ="), states[i], "\n")
    print(formalGen(product(c(m, 0)), states[i]))
    cat(noquote("EVGenerator \t\tdiscVar ="), states[i], "\n")
    print(EVGen(states[i]))
  }
  
  # sum over all discrete states:
  pg <- Reduce("+", lapply(1:k, function(i){
   lone(n+i, n+k)*increase_arity(polyGen(states[i]), n+1:k)
  }))
  evg <- Reduce("+", lapply(1:k, function(i) EVGen(states[i])))
  cat("\nsum over all discVars:\n")
  cat("polyGenerator\t")
  print(pg)
  cat("EVGenerator\t")
  print(evg)
  cat("difference\t ")
  print(pg-evg)
}
#------------- tests --------------------------

test_that("generator works for model 1", {
  
  #### definitions
  data("genePdmp1")
  data("genePoly1")
  parms(genePdmp1)["κ10"] <- 2
  parms(genePoly1)["κ10"] <- 2
  n <- length(genePoly1@init) - 1
  states <- genePoly1@discStates[[1]]
  k <- length(states)
  
  #------- compare with known generator -------
  
  f <- 5*linear(1:2)
  
  formalGen <- function(f, discVar){ # definition of known formula
    gen <- quote(
      deriv(f,1)*linear(c(-β,α)) + 
        (discVar*κ10-(1-discVar)*κ01)*funcdiff(f,0,1)
    )
    gen <- with(as.list(parms(genePoly1)), eval(gen))
    return(subs(gen, n+1, discVar))
  }
  
  expect_identical(polyGenerator(genePoly1)(f)(1), 
                   formalGen(f, 1))
  
  #------ compare with EVGenerator -----------
  
  m <- 3
  mp <- product(c(m, 0))
  #printGenerators(genePoly1, formalGen, m)
  
  # sum over all discrete states
  a <- Reduce("+", lapply(1:k, function(i){
    gen <- polyGenerator(genePoly1)(mp)(states[i])
    lone(n+i, n+k)*increase_arity(gen, n+1:k)
  }))
  b <- Reduce("+", lapply(1:k, function(i)
    EVGenerator(genePoly1, m, states[i])
  ))
  expect_equal(a, b)
  
  #------- compare with generator from pdmpsim ---------
  
  fvals <- seq(from = 0, to = 10, by = 1)
  g <- function(ξ, θ) 5*ξ + 10*θ
  gp <- 5*linear(1:2)
  
  pdmpGen <- function(discVar, fvals){
    sapply(fvals, function(val)
      generator(genePdmp1)(g)(t = 1, x = c("ξ" = val, "θ" = discVar))
    )
  }
  
  polyGen <- function(discVar, fvals){
    sapply(fvals, function(val){
      as.function.spray(formalGen(gp, discVar))(val)
    }
    )
  }
  
  test1 <- cbind(pdmpGen(0, fvals), pdmpGen(1, fvals)) 
  test2 <- cbind(polyGen(0, fvals), polyGen(1, fvals)) 
  expect_equal(test1, test2, check.attributes = FALSE)
})

test_that("generator works for model 2", {
  
  #### definitions
  data("genePdmp2")
  data("genePoly2")
  n <- length(genePoly2@init) - 1
  states <- genePoly2@discStates[[1]]
  k <- length(states)
  
  #------- compare with known generator -------
  
  f <- product(1:3)
  
  formalGen <- function(f, discVar){ # definition of known formula
    gen <- quote(
      deriv(f,1)*linear(c(-β1, 0, α1)) + deriv(f,2)*linear(c(α2, -β2, 0)) + 
        switch(discVar+1, -κ01, κ10)*funcdiff(f,0,1)
    )
    gen <- with(as.list(parms(genePoly2)), eval(gen))
    return(subs(gen, n+1, discVar))
  }

  expect_identical(polyGenerator(genePoly2)(f)(1), 
                   formalGen(f, 1))
  
  #------ compare with EVGenerator -----------
  
  m <- c(2, 10)
  mp <- product(c(m, 0))
  #printGenerators(genePoly2, formalGen, m)
  
  # sum over all discrete states
  a <- Reduce("+", lapply(1:k, function(i){
    gen <- polyGenerator(genePoly2)(mp)(states[i])
    lone(n+i, n+k)*increase_arity(gen, n+1:k)
  }))
  b <- Reduce("+", lapply(1:k, function(i) 
    EVGenerator(genePoly2, m, states[i])
  ))
  expect_equal(a, b)
  
  #------- compare with generator from pdmpsim ---------
  
  gvals <- seq(from = 0, to = 10, by = 1)
  g <- function(ξ1, ξ2, θ) 5*ξ1*ξ2^2 + θ
  gp <- 5*product(c(1,2,0)) + lone(3,3)
  
  pdmpGen <- function(discVar, gvals){
    sapply(gvals, function(val)
      generator(genePdmp2)(g)(t = 1, 
                              x = c("ξ1" = val, "ξ2" = 3*val, "θ" = discVar))
    )
  }
  
  polyGen <- function(discVar, gvals){
    sapply(gvals, function(val){
      as.function.spray(formalGen(gp, discVar))(c(val, 3*val))
    })
  }
  
  test1 <- cbind(pdmpGen(0, gvals), pdmpGen(1, gvals)) 
  test2 <- cbind(polyGen(0, gvals), polyGen(1, gvals)) 
  expect_equal(test1, test2, check.attributes = FALSE)
})

test_that("generator works for the toggleSwitch model", {
  
  #### definitions
  data("genePdmp7")
  data("genePoly7")
  n <- length(genePoly7@init) - 1
  states <- genePoly7@discStates[[1]]
  k <- length(states)
  
  #------- compare with known generator -------
  
  f <- product(1:3)
  
  # formalGen works only for sprays with funcdiff(...) != 0
  formalGen <- function(f, discVar){ 
    gen <- quote(
      deriv(f,1)*(-βA*lone(1,3) + switch(discVar, 0, αA, 0, αA)) + 
      deriv(f,2)*(-βB*lone(2,3) + switch(discVar, 0, 0, αB, αB)) +
      switch(discVar,
               κ01B*funcdiff(f,3,1)+κ01A*funcdiff(f,2,1),
               κ01B*funcdiff(f,4,2) + lone(2,3)*funcdiff(f,1,2),
               κ10B*lone(1,3)*funcdiff(f,1,3)+κ01A*funcdiff(f,4,3),
               κ10B*lone(1,3)*funcdiff(f,2,4)+κ10A*lone(2,3)*funcdiff(f,3,4)
      )
    )
    gen <- with(as.list(parms(genePoly7)), eval(gen))
    return(subs(gen, n+1, discVar))
  }
  
  expect_identical(polyGenerator(genePoly7)(f)(1), 
                   formalGen(f, 1))
  
  #------ compare with EVGenerator -----------
  
  m <- c(4, 5)
  mp <- product(c(m, 0))
  #printGenerators(genePoly7, formalGen, m)
  
  # sum over all discrete states
  a <- Reduce("+", lapply(1:k, function(i){
    gen <- polyGenerator(genePoly7)(mp)(states[i])
    lone(n+i, n+k)*increase_arity(gen, n+1:k)
  }))
  b <- Reduce("+", lapply(1:k, function(i) 
    EVGenerator(genePoly7, m, states[i])
  ))
  expect_true(a == b)
  
  #------- compare with generator from pdmpsim ---------
  
  gvals <- seq(from = 0, to = 10, by = 1)
  g <- function(ξA, ξB, θ) 5*ξA*ξB^2 + θ
  gp <- 5*product(c(1,2,0)) + lone(3,3)
  
  pdmpGen <- function(discVar, gvals){
    sapply(gvals, function(val)
      generator(genePdmp7)(g)(t = 1, 
                              x = c("ξA" = val, "ξB" = 6*val, "θ" = discVar))
    )
  }
  
  polyGen <- function(discVar, gvals){
    sapply(gvals, function(val){
      as.function.spray(formalGen(gp, discVar))(c(val, 6*val))
    })
  }
  
  test1 <- cbind(pdmpGen(1, gvals), pdmpGen(3, gvals), pdmpGen(4, gvals)) 
  test2 <- cbind(polyGen(1, gvals), polyGen(3, gvals), polyGen(4, gvals)) 
  expect_equal(test1, test2, check.attributes = FALSE)
})

test_that("generator works for constant polynomials", {
  
  #### definitions
  data("genePdmp2")
  data("genePoly2")
  n <- length(genePoly2@init) - 1
  states <- genePoly2@discStates[[1]]
  k <- length(states)
  
  m <- c(0, 0)
  mp <- product(c(m, 0))
  
  expect_true(is.zero(polyGenerator(genePoly2)(mp)(1)))
  
  b <- Reduce("+", lapply(1:k, function(i){
    EVGenerator(genePoly2, m, states[i])}
  ))
  expect_true(is.zero(b))
})

###### alt #########

test_that("old code not suited for a test", {
  skip("")

# define generator for every polyModel (from formal definition)
formalGenerator <- function(nr, f){
  
  polyMod <- get(paste("polyModel", nr, sep = '')) # polyPdmp-model
  funcdiff <- function(f, val1, val2){ # a help function
    subs(f, arity(f), val1, TRUE)-subs(f, arity(f), val2, TRUE)
  }
  
  function(discVar){
    with(as.list(polyMod@parms),
         switch(nr,
                #models 1-6:
                deriv(f,1)*linear(c(-β,α)) + (discVar*κ10-(1-discVar)*κ01)*funcdiff(f,0,1),
                deriv(f,1)*linear(c(-β1, 0, α1)) + deriv(f,2)*linear(c(α2, -β2, 0)) + switch(discVar+1, -κ01, κ10)*funcdiff(f,0,1), # former modell 9
                deriv(f,1)*linear(c(-β,α)) + switch(discVar+1, -κ01, κ10*lone(1,2))*funcdiff(f,0,1),
                deriv(f,1)*linear(c(-β,α)) + switch(discVar+1, -κ01*lone(1,2), κ10)*funcdiff(f,0,1), # former modell 2
                deriv(f,1)*(-β*lone(1,2) + switch(discVar+1, α0, α1)) + switch(discVar+1, -κ01*lone(1,2), κ10)*funcdiff(f,0,1),
                deriv(f,1)*linear(c(-β,α)) + switch(discVar+1, -κ01p*lone(1,2) - κ01c*one(2), κ10p+κ10c)*funcdiff(f,0,1),
                #model 7: #(nicht mehr aktuell)
                deriv(f,1)*(-βA*lone(1,3) + switch(discVar, 0, αA, 0, αA)) + deriv(f,2)*(-βB*lone(2,3)+switch(discVar, 0, 0, αB, αB))
                + switch(discVar, 
                         κ01A*funcdiff(f,3,1)+κ01B*funcdiff(f,2,1),
                         κ01A*funcdiff(f,4,2)+κ10B*lone(1,3)*funcdiff(f,1,2),
                         κ10A*lone(2,3)*funcdiff(f,1,3)+κ01B*funcdiff(f,4,3),
                         κ10A*lone(2,3)*funcdiff(f,2,4)+κ10B*lone(1,3)*funcdiff(f,3,4)
                ),
                #model 8:
                deriv(f,1)*(linear(c(-2*κ01, 0, 0),2) + linear(c(-β,2*κ10,α))) + deriv(f,2)*(linear(c(κ01, 0, 0),2) - κ10*lone(2,3))
                + switch(discVar+1, -κ01*lone(2,3), κ10)*funcdiff(f,0,1),
         )
    )
  }
}

})