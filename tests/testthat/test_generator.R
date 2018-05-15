#======== todo =================================================================

library(spray)
context("generator")






###### alt #########

test_that("old code not suited for a test", {
  skip()

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