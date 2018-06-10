### random examples ###

unimodal.example <- function(A=0, B=10){
  dis <- random.distribution(A, B)
  is.unimodal(A, B, mmix(1:4, A, B, distrib = dis))
}
unimodal.example()

### examples that do not work well ###

d1 <- list(list(spec = "unif", min = 3, max = 4),
           list(spec = "lnorm", meanlog = 4))
w1 <- c(10,1)
curve(dmix(a=0, b=10, w1, distrib = d1)(x), 0, 11)
is.unimodal(0, 10, mmix(1:4, 0, 10, w1, distrib = d1))

d2 <- list(list(spec = "unif", min = 6, max = 7),
           list(spec = "norm", mean = 3))
w2 <- c(1,2)
curve(dmix(a=0, b=8, w2, distrib = d2)(x), 0, 8)
is.unimodal(0, 8, mmix(1:4, 0, 8, w2, distrib = d2))

d3 <- list(list(spec = "unif", min = 0, max = 1),
           list(spec = "unif", min = 2, max = 3))
w3 <- c(1,1) # works with c(1,1.2)
B <- 4       # works with B = 3
curve(dmix(a=0, B, w3, distrib = d3)(x), -1, B+1)
is.unimodal(0, B, mmix(1:4, 0, B, w3, distrib = d3))
