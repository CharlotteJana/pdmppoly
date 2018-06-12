### random example ###

  unimodal.example <- function(A = 0, B = 10){
    dis <- random.distribution(A, B)
    is.unimodal(A, B, mmix(1:4, A, B, distrib = dis))
  }
  unimodal.example()


### mixture of two normal distributions ###

  d <- list(list("norm", mean = 3),
            list("norm", mean = 13))
  
  # recognized as not unimodal:
  is.unimodal(0, 20, mmix(a = 0, b = 20, distrib = d, weights = c(0.3, 0.7)))
  # bimodal but not recognized as such:
  is.unimodal(0, 20, mmix(a = 0, b = 20, distrib = d))
  
  
  # bimodal but not recognized as such:
  d[[2]]["mean"] <- 8
  is.unimodal(0, 20, mmix(a = 0, b = 20, distrib = d, weights = c(0.3, 0.7)))
  

### examples that do not work well ###

  # very close to unimodal:
  d1 <- list(list(spec = "unif", min = 3, max = 4),
             list(spec = "lnorm", meanlog = 4))
  w1 <- c(10,1)
  curve(dmix(a = 0, b = 10, weights = w1, distrib = d1)(x), 0, 11)
  is.unimodal(0, 10, mmix(1:4, 0, 10, weights = w1, distrib = d1)) # FALSE

  # bimodal but not recognized as such:
  d2 <- list(list(spec = "unif", min = 6, max = 7),
             list(spec = "norm", mean = 3))
  w2 <- c(1,2)
  curve(dmix(a = 0, b = 8, weights = w2, distrib = d2)(x), 0, 8)
  is.unimodal(0, 8, mmix(1:4, 0, 8, weights = w2, distrib = d2)) # TRUE

