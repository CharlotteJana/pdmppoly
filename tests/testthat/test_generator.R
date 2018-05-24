#======== todo =================================================================

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
  EVGen <- function(discVarIndex) 
    EVGenerator(polyMod, m, discVarIndex)
  
  for(i in 1:k){
    cat(noquote("polyGenerator    \tdiscVar ="), states[i], "\n")
    print(polyGen(states[i]))
    cat(noquote("formal generator \tdiscVar ="), states[i], "\n")
    print(formalGen(states[i]))
    cat(noquote("EVGenerator \t\tdiscVar ="), states[i], "\n")
    print(EVGen(i))
  }
  
  # sum over all discrete states:
  pg <- Reduce("+", lapply(1:k, function(i){
   lone(n+i, n+k)*increase_arity(polyGen(states[i]), n+1:k)
  }))
  evg <- Reduce("+", lapply(1:k, function(i) EVGen(i)))
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
    EVGenerator(genePoly1, m, i)
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
                #model 7:
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

compareGenerators <- function(nr, x, ...){
  UseMethod("compareGenerators", x)
}

compareGenerators.spray <- function(nr, f){
  
  # get polyPdmp-model
  polyMod <- get(paste("polyModel", nr, sep = ''))
  n <- length(polyMod@init)-1
  dom <- polyMod@discDomain
  
  # get polyGenerator and formalGenerator
  polyGen <- polyGenerator(polyMod)(f) # function(discVar)
  formalGen <- formalGenerator(nr, f) # function(discVar)
  
  # compare them
  for(i in 1:length(dom)){
    cat(noquote("polyGenerator    \tdiscVar ="), dom[i], "\n")
    print(polyGen(dom[i]))
    cat(noquote("formal generator \tdiscVar ="), dom[i], "\n")
    print(subs(formalGen(dom[i]), n+1, dom[i]))
  }
}

compareGenerators.numeric <- function(nr, m){
  
  # definitions
  f <- product(c(m, 0)) # f := z₁ᵐ¹*...*zₙᵐⁿ with arity = n+1
  polyMod <- get(paste("polyModel", nr, sep = ''))
  n <- length(polyMod@init) - 1
  dom <- polyMod@discDomain
  polyGen <- polyGenerator(polyMod)(f) # function(discVar)
  formalGen <- formalGenerator(nr, f) # function(discVar)
  
  if(length(m) != n) stop("m has to be a vector of length ", n)
  
  # compare them
  for(i in 1:length(dom)){
    cat(noquote("polyGenerator    \tdiscVar ="), dom[i], "\n")
    print(blowupSpray(polyMod, polyGen(dom[i])))
    cat(noquote("formal generator \tdiscVar ="), dom[i], "\n")
    print(blowupSpray(polyMod, subs(formalGen(dom[i]), n+1, dom[i])))
    cat(noquote("EVGenerator \t\tdiscVar ="), dom[i], "\n")
    print(EVGenerator(polyMod, m, i))
  }
  
  # sum over all discrete variables:
  pg <- Reduce("+", lapply(1:length(dom), function(i){
    lone(n+i, n+length(dom))*blowupSpray(polyMod, polyGen(dom[i]))
  }))
  evg <- Reduce("+", lapply(1:length(dom), function(i) EVGenerator(polyMod, m, i)))
  cat("\nsum over all discVars:\n")
  cat("polyGenerator\t")
  print(pg)
  cat("EVGenerator\t")
  print(evg)
  cat("difference\t ")
  print(pg-evg)
  
  invisible()
}
})