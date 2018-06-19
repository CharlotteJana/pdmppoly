#======== todo =================================================================

# # nur in dieser datei benutzt, sonst nicht
# compute.central.moment <- function(j, moments){
#   e <- moments[1]
#   s <- sapply(1:j, function(k) choose(j,k)*moments[k]*(-e)^(j-k))
#   print(s)
#   sum(s) + (-e)^j
# }
# 
# # nirgends benutzt
# centr.moments <- function(moments){
#   sapply(1:length(moments), function(j) compute.central.moment(j, moments))
# }
# 
# # nirgends benutzt, vllt = standardize.moments?
# stand.moments <- function(moments){
#   sd <- sqrt(moments[2]-moments[1]^2)
#   sapply(1:length(moments), function(j) compute.central.moment(j,moments)/sd^j)
# }
# 

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