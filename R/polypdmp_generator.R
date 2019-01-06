#======== todo =================================================================
#t2 warum ist generator eine funktion von discvar????
#t2 EVGenerator umbenennen
#t2 EVGenerator und Generator unabhängig von der Stelle von θ in init 
#t3 arity sollte eigentlich n sein?

#' @include polypdmp_class.R polypdmp_accessors.R
NULL

#' Generator
#' 
#' @inherit pdmpsim::generator description
#' @param obj an object of class \code{\link{polyPdmpModel}}.
#' @return The generator \code{Q} of \code{obj} as defined above. This is a
#'   function which takes as argument a single polynomial \code{f} (represented
#'   as spray object). The variables of \code{f} represent the variables of the
#'   PDMP, given in the same order as the variables given in slot \code{init(obj)}. 
#'   The result  
#'   \code{Q(f)} is a function with parameter \code{discVar}. As can be seen
#'   in the formula, the resulting polynomial \code{Q(f)(i, z)} depends on the value
#'   \code{i} of the discrete variable. The returned value of function \code{Q(f)}
#'   is a polynomial represented as \code{spray} object.
#' @examples
#' library(spray)
#' data("simplePoly")
#' g <- product(c(1,1))
#' generator(simplePoly)(g)(-1)
#' 
#' # comparison with theoretic solution:
#' Qg_theoretic <- product(c(2,0))-2*product(c(1,1))
#' identical(generator(simplePoly)(g)(1), subs(Qg_theoretic, 1, 1))
#' @note Method \code{generator} only works for one discrete variable and
#' this variable should be the last entry in slot \code{init}.
#' @importFrom spray arity deriv subs lone 
#' @export
setMethod("generator", signature(obj = "polyPdmpModel"), function(obj) {
 function(poly){
   function(discVar){
     
     ### definitions
     # discVar = value of the discrete variable
     n <- length(obj@init) - length(obj@discStates) # number of continuous variables
     nj <- length(obj@ratesprays) # number of jumptypes
     stateIndex <- getIndex(discVar, obj@discStates[[1]])
     
     ### poly
     if(arity(poly) != n+1) stop("arity of the polynomial has to equal ", n+1) # arity sollte eigentlich n sein?
     
     ### sum over continuous states
     list <- lapply(1:n, function(i){
       ifelse(is.zero(deriv(poly, i)),
              return(0*lone(1, arity(poly))),
              return(deriv(poly,i)*obj@dynsprays[[i]][[stateIndex]]))
       })
     s1 <- Reduce("+", list)
     
     # sum over possible new discrete states
     ratematrix <- ratespraysToMatrix(obj)
     list <- lapply(seq_along(obj@discStates[[1]]), function(j){
       if(!is.null(ratematrix[[stateIndex]][[j]])){
         diff <- increase_arity(subs(poly, n+1, obj@discStates[[1]][j]), n+1) - 
           increase_arity(subs(poly, n+1, discVar), n+1)
         if(!is.zero(diff))
           return(ratematrix[[stateIndex]][[j]] * diff)
        }
       return(0*lone(1, n+1))
       })
     s2 <- Reduce("+", list)
     
     return(subs(s1+s2, n+1, discVar)) # increase_arity???
   }
 } 
})

#' Generator of m-th moment
#' 
#' Let Q be the \code{\link{generator}} of a PDMP. Method \code{EVGenerator}
#' computes the generator Q(f) where f is given as \deqn{f(i,z) = δ_{ji}\cdot
#' z^{m[1]}z^{m[2]} ... z^{m[k]}}{f(i,z) = δᵢⱼ*z₁ᵐ¹*...*zₙᵐⁿ } with a given
#' vector \code{m} and \code{j} being a fixed state of the discrete variable d.
#' This method is needed to compute the expected value \eqn{\mathbb{E}(X_t^m | d
#' = j)}{E(Xₜᵐ | d=j)}.
#' @param obj object of class \code{\link{polyPdmpModel}}.
#' @param j integer giving the state of the discrete variable.
#' @param m integer vector whose length equals the number of continous
#'   variables.
#' @return a polynomial represented as \code{spray} object. This polynomial has
#'   more variables than \code{init(obj)} because the discrete variable \code{d}
#'   is replaced by different indicator variables \code{d1, ..., dk}, one for
#'   every possible state. See \code{\link{blowupSpray}} for more details.
#' @seealso \code{\link{generator}} to compute the generator and obtain the
#'   result as spray object with the same number of variables as
#'   \code{init(obj)}, \code{\link{generator}} to compute the generator and
#'   obtain the result as function.
#' @note This method only works for one discrete variable. This variable has to
#' be the last element of vector \code{init} of model \code{obj}.
#' @importFrom spray product deriv subs lone
EVGenerator <- function(obj, m, j){ 
  # Computes EV(generator(θⱼ*spray)), where
  # θⱼ is the i-th indicator variable for the discrete variable θ (j ϵ {1,...,k}),
  # spray = z₁ᵐ¹*...*zₙᵐⁿ with z = vector of all continous variables.
  # The result is a (blown up) spray of arity n+k.
  
  ### definitions
  index <- getIndex(j, obj@discStates[[1]])
  k <- length(obj@discStates[[1]]) # number of different discrete states
  n <- length(obj@init) - 1 # number of continuous variables
  if(length(m) != n) stop("length(m) has to equal ", n)
  spray <- product(c(m, 0)) # spray := z₁ᵐ¹*...*zₙᵐⁿ with arity = n+1
  ratematrix <- ratespraysToMatrix(obj)
  
  ### sum over continuous states r: ∂spray/∂zᵣⱼ * dynfunc(r) * θⱼ
  list <- lapply(1:n, function(r){
    ifelse(is.zero(deriv(spray, r)),
           return(0*lone(1, arity(spray))),
           return(deriv(spray,r)*obj@dynsprays[[r]][[index]]))
  })
  s1 <- Reduce("+", list) # sum
  s1 <- increase_arity(subs(s1, n+1, j), n+1:k) # substitute θ with j
  if(!is.zero(s1))
    s1 <- s1*lone(n+index, n+k) # times θⱼ
  
  ### sum over q ϵ discStates: rate(q→j)*spray*θⱼ
  list <- lapply(1:k, function(q){ # rate(q→j)*θⱼ
    if(!is.null(ratematrix[[q]][[index]])) {
      h <- subs(ratematrix[[q]][[index]], n+1, obj@discStates[[1]][q])
      increase_arity(h, n+1:k)*lone(n+q,n+k)
    }
    else 0*lone(1,n+k)
  })
  s2 <- Reduce("+", list) # sum
  s2 <- s2*increase_arity(spray, n + 2:k)# times spray


  ### sum over q ϵ discStates: -rate(j→q)*spray*θⱼ
  list <- lapply(1:k, function(q){ # rate(j→q) 
    if(!is.null(ratematrix[[index]][[q]])) 
      increase_arity(subs(ratematrix[[index]][[q]], n+1, j), n+1:k)
    else 0*lone(1,n+k)
  })
  s3 <- Reduce("+", list) # sum
  s3 <- -s3*increase_arity(spray, n+2:k)*lone(n+index, n+k) # times -spray*θⱼ
  
  sum <- s1+s2+s3
  return(sum)
}


