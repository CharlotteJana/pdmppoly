#======== todo =================================================================
#t1 references for inequalities
#t1 munif?
#t2 tests
#t1 lower bound for 2-b-unimodal?
#t1 documentation of parameters

#' Test if moments come from a unimodal distribution with compact support
#'
#' Given a compact support \code{[A, B]} and moments \code{m1, m2, ...} of a 
#' one dimensional distribution, this method checks if the (unknown) 
#' distribution can be unimodal. \cr 
#' There exist several inequalities that test for nonunimodality in the
#' literature. Depending on the inequality, moments up to order 2 or 4 are
#' required. A distribution that satisfies all inequalities that contain only
#' moments up to order 2 is called \emph{2-b-unimodal}. A distribution that
#' satisfies all inequalities that contain only moments up to order 4 is called
#' \emph{4-b-unimodal}. The internal methods \code{is.2_b_unimodal} and
#' \code{is.4_b_unimodal} test these inequalities. Method \code{is.unimodal}
#' performs these checks, depending on the number of given distributions. A
#' 4-b-unimodal distribution can still be nonunimodal (see the examples below).
#' But a failure of these tests assures the distribution not to be unimodal.
#'
#' @inheritParams exists.distribution
#' @param eps numeric value. Some inequalities are of the form \code{... > 0}.
#'   For numerical reasons it is better to test for \code{... > eps} where
#'   \code{eps} is a small number.
#' @return boolean value indicating if the corresponding distribution
#' is 2-/4-b-unimodal (TRUE) or not (FALSE).
#' @example inst/examples/unimodal.R
#' @name is.unimodal
#' @export
is.unimodal <- function(lower, upper, moments, eps = 1e-10){
  #[A,B] = Support der ZG
  #m = sortierter Vektor mit Momenten
  
  if(length(moments) < 4){
    is.2_b_unimodal(lower, upper, moments, eps)
  }
  if(length(m) < 6){
    is.2_b_unimodal(lower, upper, moments, eps)
    is.4_b_unimodal(lower, upper, moments, eps)
  }
}

#' @rdname is.unimodal
#' @export
is.2_b_unimodal <- function(lower, upper, moments, eps = 1e-10){
  
  if(is.vector(moments))
    moments <- t(moments)
  
  # [a, b] = support of the standardized distribution
  sd <- sqrt(moments[, 2]-moments[, 1]^2) # standard derivation
  a <- (lower-moments[, 1])/sd 
  b <- (upper-moments[, 1])/sd
  
  # Existence of a distribution
  eqn1 <- ifelse(1+a*b <= 0, TRUE, FALSE)
  
  # Lower bound ???
  # ...
  
  # Upper bound (Teuscher, Guiard)
  eqn2 <- ifelse(a + b > eps & b > sqrt(3) & a <= -b+sqrt(b^2-3), TRUE, FALSE)
  eqn3 <- ifelse(a+b < -eps & a < -sqrt(3) & b >= -a-sqrt(a^2-3), TRUE, FALSE)
  eqn4 <- ifelse(abs(a+b) <= eps & b >= sqrt(3), TRUE, FALSE)
  eqn5 <- (eqn2 | eqn3 | eqn4)
  
  results <- dplyr::case_when(
    isTRUE(eqn5) ~ "2-b-unimodal",
    !isTRUE(eqn5) ~ "not unimodal",
    !isTRUE(eqn1) ~ "not existant",
    TRUE ~ rep(NA_character_, length(eqn1))
  )
  
  return(results)
}

#' @rdname is.unimodal
#' @export
is.4_b_unimodal <- function(lower, upper, moments, eps = 1e-10){

  # [a, b] = support of the standardized distribution
  sd <- sqrt(moments[[2]]-moments[[1]]^2) # strandard derivation
  a <- (lower-moments[[1]])/sd
  b <- (upper-moments[[1]])/sd
  g1 <- (moments[[3]]-3*sd^2*moments[[1]]-moments[[1]]^3)/sd^3
  g2 <- (-3*moments[[1]]^4+6*moments[[2]]*moments[[1]]^2-
           4*moments[[3]]*moments[[1]]+moments[[4]])/sd^4-3

  
  Q <- 4*g1*(a+b)+(3-a^2)*(3-b^2)

  if (1+a*b > 0 | g2 < g1^2-2 | g2 > g1^2-2-(a*b*(g1-a+1/a)*(g1-b+1/b))/(1+a*b) | g1 > b-1/b | g1 < a-1/a){
    message("There is no distribution with these parameters, because ")
    if(1+a*b > 0) message("1+a*b > 0.")
    if(g2 < g1^2-2) message("γ2 < γ1²-2.")
    if(g2 > g1^2-2-(a*b*(g1-a+1/a)*(g1-b+1/b))/(1+a*b)) 
      message("γ2 > γ1²-2-(a*b*(γ1-a+1/a)*(γ1-b+1/b))/(1+a*b) .")
    if(g1 > b-1/b) message("γ1 > b-1/b.")
    if(g1 < a-1/a) message("γ1 < a-1/a.")
    return(FALSE)
  }
  
  ### Upper bound (Teuscher, Guiard)
  
  if (a + b > eps){ # Fall a < b
    if (Q >= 0){
      eq <- (4*g1*(a^2+b^2+a*b-3) 
             -2*Q*(3+a*b+sqrt(Q))/(a+b))/(5*(a+b))-6/5
      boolUpper <- (g2 <= eps + eq)
    } else {
      eq <- (4*g1*(a^2+b^2+a*b-3) +
             Q*(4*g1+(a+b)*(a^2-3))/(3+2*a*b+a^2))/(5*(a+b))-6/5
      boolUpper <- (g2 <= eps + eq)
    }
  } else if (a + b < -eps){ # Fall a > b
    if (Q >= 0){
      eq <- (4*g1*(a^2+b^2+a*b-3)
             -2*Q*(3+a*b-sqrt(Q))/(a+b))/(5*(a+b))-6/5
      boolUpper <- (g2 <= eps + eq)
    } else {
      eq <- (4*g1*(a^2+b^2+a*b-3)
             +Q*(4*g1+(a+b)*(b^2-3))/(3+2*a*b+b^2))/(5*(a+b))-6/5
      boolUpper <- (g2 <= eps + eq)
    }
  } else{  # Fall a = -b
      eq <- 1/5*((12*g1^2)/(3-a^2)+3*a^2-15)
      boolUpper <- (g2 <= eps + eq)
  }

  ### Lower bound (Johnsson, Rogers)
  
  s <- sqrt(as.complex(g2+6))
  t <- 9*sqrt(as.complex(5*g2+6))*sqrt(as.complex(16*g2+21)) + s*(40*g2+51)
  calc_q <- function(z){
    ((5^(1/3)*s*t^(1/3) + Conj(z)*t^(2/3) + z*4*5^(2/3)*g2 + 
        z*6*5^(2/3))/(9*5^(1/3)*s*t^(1/3)))
  }
  # all 3 solutions of ...
  q <- c(calc_q(-1-sqrt(3)*1i), calc_q(-1+sqrt(3)*1i), calc_q(2)) 
  
  # only real solutions between 0 and 1:
  q <- Re(q[Im(q) < eps & 0 <= Re(q) & Re(q) <= 1]) 
  
  if (length(q) == 1){
    boolLower <- (g1^2 <= eps + (108*q[1]^4)/((1-q[1])*(1+3*q[1])^3))
  } else if (length(q) == 2){
      q <- sort(q)
      eq1 <- (108*q[1]^4)/((1-q[1])*(1+3*q[1])^3)
      eq2 <- (g1^2 <= eps + (108*q[2]^4)/((1-q[2])*(1+3*q[2])^3))
      boolLower <- (g1^2 + eps >= eq1 & eq2)
  } else {
      message("Error during calculation of the lower bound.")
  }
  
#   if (g2 < -1.3125){
#     print("Fall g2 < -1.3125")
#     boolLower <- FALSE
#   } else if (g2 > -1.2){ #in diesem Fall ist q reell und positiv!  Lieber >= statt >? Gilt auch q <= 1?
#       boolLower <- (g1^2 <= (108*q[1]^4)/((1-q[1])*(1+3*q[1])^3))
#       print("Fall g2 > -1.2")
#       print(length(q))
#   } else {
#     print("Fall  -1.3125 < g2 < -1.2")
#     print(length(q))
#     #f <- function(q){(72*q^2*(3*q-1))/((1-q)*(1+3*q)^2)-5*g2+6}
#     #q <- uniroot(f, c(-1,0.2), tol = eps)$root
#     
#     boolLower <- TRUE
#   }

  ### return
  
  if(boolUpper & boolLower){
    message("The distribution is 4-b-unimodal.")
  } else {    
    message("The distribution cannot be unimodal.")
  }
  return(boolUpper & boolLower)
}

##################################################################
########  Beispiele, bei denen es nicht gut funktioniert #########
##################################################################

d1 <- list(list(spec = "unif", min = 3, max = 4),
           list(spec = "lnorm", meanlog = 4))
w1 <- c(10,1)
#curve(dmix(a=0, b=10, w1, distrib = d1)(x), 0, 11)
#is.unimodal(0, 10, mmix(1:4, 0, 10, w1, distrib = d1))

d2 <- list(list(spec = "unif", min = 6, max = 7),
           list(spec = "norm", mean = 3))
w2 <- c(1,2)
#curve(dmix(a=0, b=8, w2, distrib = d2)(x), 0, 8)
#is.unimodal(0, 8, mmix(1:4, 0, 8, w2, distrib = d2))

d3 <- list(list(spec = "unif", min = 0, max = 1),
           list(spec = "unif", min = 2, max = 3))
w3 <- c(1,1) # mit c(1,1.2) gehts
B <- 4       # mit B = 3 gehts
#curve(dmix(a=0, B, w3, distrib = d3)(x), -1, B+1)
#is.unimodal(0, B, mmix(1:4, 0, B, w3, distrib = d3))
