##### generator #####

setGeneric("polyGenerator", function(obj, ...) standardGeneric("polyGenerator"))
setMethod("polyGenerator", signature(obj = "polyPdmpModel"), function(obj) {
 function(poly){
   function(discVar){
     
     ### definitions
     # discVar = value of the discrete variable
     n <- length(obj@init) - 1 # number of continuous variables
     nj <- length(obj@ratesprays) # number of jumptypes
     discDomainIndex <- getIndex(discVar, obj@discDomain)
     
     ### poly
     if(arity(poly) != n+1) stop("arity of the polynomial has to equal ", n+1) # arity sollte eigentlich n sein?
     
     ### sum over continuous states
     list <- lapply(1:n, function(i) deriv(poly,i)*obj@dynsprays[[i]][[discDomainIndex]])
     s1 <- Reduce("+", list)
     
     # sum over possible new discrete states
     ratematrix <- ratespraysToMatrix(obj)
     list <- lapply(1:length(obj@discDomain), function(j){
       if(!is.null(ratematrix[[discDomainIndex]][[j]])) ratematrix[[discDomainIndex]][[j]]*(subs(poly, n+1, obj@discDomain[j], keepArity = TRUE) - subs(poly, n+1, discVar, keepArity = TRUE))
       else 0*lone(1,n+1)
       })
     s2 <- Reduce("+", list)
     
     return(subs(s1+s2, n+1, discVar))
   }
 } 
})

EVGenerator <- function(obj, m, i){ 
  # Computes EV(generator(θᵢ*spray)), where
  # θᵢ is the i-th indicator variable for the discrete variable θ (i ϵ {1,...,k}),
  # spray = z₁ᵐ¹*...*zₙᵐⁿ with z = vector of all continous variables.
  # The result is a (blown up) spray of arity n+k.
  
  ### definitions
  dom <- obj@discDomain
  k <- length(dom) # number of different discrete states
  n <- length(obj@init) - 1 # number of continuous variables
  if(length(m) != n) stop("length(m) has to equal ", n)
  spray <- product(c(m, 0)) # spray := z₁ᵐ¹*...*zₙᵐⁿ with arity = n+1
  ratematrix <- ratespraysToMatrix(obj)
  
  ### sum over continuous states j: ∂spray/∂zⱼ * dynfunc(j) * θᵢ
  list <- lapply(1:n, function(j) deriv(spray,j)*obj@dynsprays[[j]][[i]])
  s1 <- Reduce("+", list) # sum
  s1 <- increase_arity(subs(s1, n+1, dom[i]), n+k) # substitute θ with i
  s1 <- s1*lone(n+i, n+k) # times θᵢ
  
  ### sum over j ϵ discDomain: rate(j→i)*spray*θⱼ
  list <- lapply(1:k, function(j){ # rate(j→i)*θⱼ
    if(!is.null(ratematrix[[j]][[i]])) {
      increase_arity(subs(ratematrix[[j]][[i]], n+1, dom[j]), n+k)*lone(n+j,n+k)
    }
    else 0*lone(1,n+k)
  })
  s2 <- Reduce("+", list) # sum
  s2 <- s2*increase_arity(spray, n+k)# times spray
  
  ### sum over j ϵ discDomain: -rate(i→j)*spray*θᵢ
  list <- lapply(1:k, function(j){ # rate(i→j)  
    if(!is.null(ratematrix[[i]][[j]])) increase_arity(subs(ratematrix[[i]][[j]], n+1, dom[i]), n+k)
    else 0*lone(1,n+k)
  })
  s3 <- Reduce("+", list) # sum
  s3 <- -s3*increase_arity(spray, n+k)*lone(n+i, n+k) # times -spray*θᵢ
  
  sum <- s1+s2+s3
  return(sum)
}  # EV(generator(θᵢ*z₁ᵐ¹*...*zₙᵐⁿ))


