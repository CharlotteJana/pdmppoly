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
    generator(polyMod)(product(c(m, 0)))(discVar)
  EVGen <- function(discVar) 
    EVGenerator(polyMod, m, discVar)
  
  for(i in 1:k){
    cat(noquote("generator    \tdiscVar ="), states[i], "\n")
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
  cat("generator\t")
  print(pg)
  cat("EVGenerator\t")
  print(evg)
  cat("difference\t ")
  print(pg-evg)
}
#------------- tests --------------------------

test_that("generator works for model K", {
  skip_if_not_installed("Deriv")
  
  #### definitions
  data("genePdmpK")
  data("genePolyK")
  parms(genePdmpK)["k10"] <- 2
  parms(genePolyK)["k10"] <- 2
  n <- length(genePolyK@init) - 1
  states <- genePolyK@discStates[[1]]
  k <- length(states)
  
  #------- compare with known generator -------
  
  f <- 5*linear(1:2)
  
  formalGen <- function(f, discVar){ # definition of known formula
    gen <- quote(
      deriv(f,1)*linear(c(-b,a)) + 
        (discVar*k10-(1-discVar)*k01)*funcdiff(f,0,1)
    )
    gen <- with(as.list(parms(genePolyK)), eval(gen))
    return(subs(gen, n+1, discVar))
  }
  
  expect_identical(generator(genePolyK)(f)(1), 
                   formalGen(f, 1))
  
  #------ compare with EVGenerator -----------
  
  m <- 3
  mp <- product(c(m, 0))
  #printGenerators(genePolyK, formalGen, m)
  
  # sum over all discrete states
  a <- Reduce("+", lapply(1:k, function(i){
    gen <- generator(genePolyK)(mp)(states[i])
    lone(n+i, n+k)*increase_arity(gen, n+1:k)
  }))
  b <- Reduce("+", lapply(1:k, function(i)
    EVGenerator(genePolyK, m, states[i])
  ))
  expect_equal(a, b)
  
  #------- compare with generator from pdmpsim ---------
  
  fvals <- seq(from = 0, to = 10, by = 1)
  g <- function(f, d) 5*f + 10*d
  gp <- 5*linear(1:2)
  
  pdmpGen <- function(discVar, fvals){
    sapply(fvals, function(val)
      generator(genePdmpK)(g)(t = 1, x = c("f" = val, "d" = discVar))
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

test_that("generator works for model K2", {
  skip_if_not_installed("Deriv")
  
  #### definitions
  data("genePdmpK2")
  data("genePolyK2")
  n <- length(genePolyK2@init) - 1
  states <- genePolyK2@discStates[[1]]
  k <- length(states)
  
  #------- compare with known generator -------
  
  f <- product(1:3)
  
  formalGen <- function(f, discVar){ # definition of known formula
    gen <- quote(
      deriv(f,1)*linear(c(-b1, 0, a1)) + deriv(f,2)*linear(c(a2, -b2, 0)) + 
        switch(discVar+1, -k01, k10)*funcdiff(f,0,1)
    )
    gen <- with(as.list(parms(genePolyK2)), eval(gen))
    return(subs(gen, n+1, discVar))
  }

  expect_identical(generator(genePolyK2)(f)(1), 
                   formalGen(f, 1))
  
  #------ compare with EVGenerator -----------
  
  m <- c(2, 10)
  mp <- product(c(m, 0))
  #printGenerators(genePolyK2, formalGen, m)
  
  # sum over all discrete states
  a <- Reduce("+", lapply(1:k, function(i){
    gen <- generator(genePolyK2)(mp)(states[i])
    lone(n+i, n+k)*increase_arity(gen, n+1:k)
  }))
  b <- Reduce("+", lapply(1:k, function(i) 
    EVGenerator(genePolyK2, m, states[i])
  ))
  expect_equal(a, b)
  
  #------- compare with generator from pdmpsim ---------
  
  gvals <- seq(from = 0, to = 10, by = 1)
  g <- function(f1, f2, d) 5*f1*f2^2 + d
  gp <- 5*product(c(1,2,0)) + lone(3,3)
  
  pdmpGen <- function(discVar, gvals){
    sapply(gvals, function(val)
      generator(genePdmpK2)(g)(t = 1, 
                              x = c("f1" = val, "f2" = 3*val, "d" = discVar))
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
  skip_if_not_installed("Deriv")
  
  #### definitions
  data("genePdmpT")
  data("genePolyT")
  n <- length(genePolyT@init) - 1
  states <- genePolyT@discStates[[1]]
  k <- length(states)
  
  #------- compare with known generator -------
  
  f <- product(1:3)
  
  # formalGen works only for sprays with funcdiff(...) != 0
  formalGen <- function(f, discVar){ 
    gen <- quote(
      deriv(f,1)*(-bA*lone(1,3) + switch(discVar, 0, aA, 0, aA)) + 
      deriv(f,2)*(-bB*lone(2,3) + switch(discVar, 0, 0, aB, aB)) +
      switch(discVar,
               k01B*funcdiff(f,3,1)+k01A*funcdiff(f,2,1),
               k01B*funcdiff(f,4,2) + lone(2,3)*funcdiff(f,1,2),
               k10B*lone(1,3)*funcdiff(f,1,3)+k01A*funcdiff(f,4,3),
               k10B*lone(1,3)*funcdiff(f,2,4)+k10A*lone(2,3)*funcdiff(f,3,4)
      )
    )
    gen <- with(as.list(parms(genePolyT)), eval(gen))
    return(subs(gen, n+1, discVar))
  }
  
  expect_identical(generator(genePolyT)(f)(1), 
                   formalGen(f, 1))
  
  #------ compare with EVGenerator -----------
  
  m <- c(4, 5)
  mp <- product(c(m, 0))
  #printGenerators(genePolyT, formalGen, m)
  
  # sum over all discrete states
  a <- Reduce("+", lapply(1:k, function(i){
    gen <- generator(genePolyT)(mp)(states[i])
    lone(n+i, n+k)*increase_arity(gen, n+1:k)
  }))
  b <- Reduce("+", lapply(1:k, function(i) 
    EVGenerator(genePolyT, m, states[i])
  ))
  expect_true(a == b)
  
  #------- compare with generator from pdmpsim ---------
  
  gvals <- seq(from = 0, to = 10, by = 1)
  g <- function(fA, fB, d) 5*fA*fB^2 + d
  gp <- 5*product(c(1,2,0)) + lone(3,3)
  
  pdmpGen <- function(discVar, gvals){
    sapply(gvals, function(val)
      generator(genePdmpT)(g)(t = 1, 
                              x = c("fA" = val, "fB" = 6*val, "d" = discVar))
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
  data("genePdmpK2")
  data("genePolyK2")
  n <- length(genePolyK2@init) - 1
  states <- genePolyK2@discStates[[1]]
  k <- length(states)
  
  m <- c(0, 0)
  mp <- product(c(m, 0))
  
  expect_true(is.zero(generator(genePolyK2)(mp)(1)))
  
  b <- Reduce("+", lapply(1:k, function(i){
    EVGenerator(genePolyK2, m, states[i])}
  ))
  expect_true(is.zero(b))
})