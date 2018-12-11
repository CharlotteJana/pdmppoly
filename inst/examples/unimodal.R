### random example ###

  unimodal.example <- function(lower = 0, upper = 10){
    dis <- random.distribution(lower, upper)
    is.unimodal(lower, upper, mmix(1:4, lower, upper, distrib = dis))
  }
  unimodal.example()

### mixture of two normal distributions ###

  d <- list(list("norm", mean = 3),
            list("norm", mean = 13))
  
  is.unimodal(0, 20, mmix(1:4, 0, 20, distrib = d, weights = c(0.3, 0.7)))

  d[[2]]["mean"] <- 8
  is.unimodal(0, 20, mmix(1:4, 0, 20, distrib = d, weights = c(0.3, 0.7)))
  

### very close to unimodal ###

  d <- list(list(spec = "unif", min = 3, max = 4),
             list(spec = "lnorm", meanlog = 4))
  w <- c(10,1)
  curve(dmix(lower = 0, upper = 10, weights = w, distrib = d)(x), 0, 11)
  is.unimodal(0, 10, mmix(1:4, 0, 10, weights = w, distrib = d))
  
  
### not unimodal at all ###
  
  d <- list(list(spec = "exp"),
            list(spec = "unif", min = 0, max = 1),
            list(spec = "unif", min = 5, max = 6),
            list(spec = "unif", min = 0, max = 10))
  curve(dmix(lower = 0, upper = 10, distrib = d)(x), -1, 11)
  is.unimodal(0,10, mmix(1:4, 0, 10, distrib = d))

