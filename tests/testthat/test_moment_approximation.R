#======== todo =================================================================
#t2 tests mit verschiedenen closure methoden

context("moment approximation")

test_that("moment calculation leads to same results for model K and model K2", {
  data(genePolyK)
  data(genePolyK2)
  init(genePolyK) <- c("ξ" = 0.5, "θ" = 1)
  init(genePolyK2) <- c("ξ1" = 0.5, "ξ2" = 0, "θ" = 1)
  parms(genePolyK) <- list(β = 0.5, α = 1, κ10 = 0.1, κ01 = 0.3)
  parms(genePolyK2) <- list(β1 = 0.5, α1 = 1, κ10 = 0.1, κ01 = 0.3, β2 = 1, α2 = 0.5)
  res1 <- momApp(genePolyK)
  res2 <- momApp(genePolyK2)
  
  expect_equal(res1$discRes, res2$discRes)
  expect_equal(res1$contRes[, 1:5], res2$contRes[, 1:5], 
               check.names = FALSE, check.attributes = FALSE)
  expect_equal(res1$moments, res2$moments[, -4], check.names = FALSE)
  
})

test_that("elements contRes, discRes and moments contain the same results", {
  data(genePolyK2)
  x <- momApp(genePolyK2, maxOrder = 6, closure = "reduceDegree")
  
  # calculate moments out of contRes and discRes
  l <- x$maxOrder
  k <- ncol(x$discRes) - 1
  n <- length(x$model@init) - length(x$model@discStates)
  names <- c("ξ1", "ξ2")
  dname <- "θ"
  
  r <- data.frame(time = rep(100, l), order = 1:l)
  
  for(i in 1:n){ # fill the moments of continuous variables
    r[, names[i]] <- rep(0, l)
    for(j in 1:l){
      ind <- rep(0, n)
      ind[i] <- j
      r[j, names[i]] <- x$contRes[nrow(x$contRes),
                                  prodlim::row.match(ind, as.matrix(x$contInd))]
    }
  }
  for(j in 1:l){ # fill the moments of the discrete variable
    r[j, dname] <- as.vector(x$model@discStates[[1]]^j %*% x$discRes[nrow(x$discRes), 2:ncol(x$discRes)])
  }
  
  expect_equal(r, x$moments[which(x$moments$time == 100), ],
               check.attributes = FALSE)
})

test_that("momApp works for a model with more than 2 discrete states", {
  data(simplePoly)
  res <- momApp(simplePoly, maxOrder = 1)
  zeros <- rep(0, nrow(res$moments))
  expect_equal(res$moments[["f"]], zeros)
  expect_equal(res$moments[["d"]], zeros)
})

test_that("order of variables in init doesn't matter", {
  skip("functionality not implemented yet")
  data(genePolyF)
  res1 <- momApp(genePolyF, closure = "reduceDegree")
  
  model <- new("polyPdmpModel",
                   descr = "Model F with different order of variables",
                   parms = parms(genePolyF), 
                   init = rev(init(genePolyF)), 
                   discStates = list(θ = 0:1),
                   dynpolys = quote(list(
                     list(overall = linear(c(-β,α)))
                   )),
                   ratepolys = quote(list(  
                     list(κ01*lone(2,2), κ10)
                   )),
                   jumpfunc = function(t, x, parms, jtype) {
                     c(1 - x[1], x[2])
                   }, 
                   times = times(genePolyF), 
                   solver = "lsodar")
  res2 <- momApp(model, closure = "reduceDegree")
  expect_identical(res1$moments, res2$moments)
  expect_identical(res1$discRes, res2$discRes)
  expect_identical(res1$contRes, res2$contRes)
})

test_that("moment calculation works for model 1", {
  
  ### definitions
  data(genePolyK)
  times(genePolyK) <- c(from = 0, to = 2000, by = 1)
  states <- discStates(genePolyK)[[1]]
  l <- 4 # only works for l < 5. If l = 5, set times["to"] = 3000
  k <- length(states)
  n <- length(genePolyK@init) - 1
  
  ### moment approximation
  momApp <- momApp(genePolyK, maxOrder = l)
  last <- nrow(momApp$contRes)
  
  ### create matrix with all moment combinations we are interested in
  s <- NULL
  for(i in 1:n){ 
    m <-  matrix(data = 0, nrow = l, ncol = n+3)
    m[,i] <- 1:l
    s <- rbind(s, m)
  }
  s <- rbind(s, c(rep(0, n), 1))
  colnames(s) <- c(names(genePolyK@init), "calc", "approx")
  s <- as.data.frame(s)
  
  ### continous variables
  for(i in 1:(n*l)){
    m <- s[i,1]
    # theoretical values:
    s[i, n+2] <- with(as.list(genePolyK@parms),
      (α/β)^m*Reduce("*", sapply(0:(m-1), function(j) (κ01+j*β)/(κ01+κ10+j*β))))
    # calculated values:
    s[i, n+3] <- momApp$contRes[last, 
                                prodlim::row.match(as.vector(s[i,1:n]), 
                                                   as.matrix(momApp$contInd))]
  }
  
  ### discrete variables
  # theoretical values:
  s[n*l+1, n+2] <- with(as.list(genePolyK@parms), κ01/(κ01+κ10))
  # calculated values:
  s[n*l+1, n+3] <- states %*% momApp$discRes[last, 2:ncol(momApp$discRes)]
 
  #print(s, digits = 4)
  for(i in 1:nrow(s)){
    expect_equal(s[i, 3], s[i, 4], tolerance = 1e-04)
  }
})

