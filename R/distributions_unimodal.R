#======== todo =================================================================
#s1 references for inequalities
#t2 NaNs zulassen

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
  if(length(moments) >= 4){
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

  # Unimodality (Teuscher, Guiard)
  eqn2 <- ifelse(a + b > eps  & b > sqrt(3),   a <= -b+sqrt(b^2-3), FALSE)
  eqn3 <- ifelse(a + b < -eps & a < -sqrt(3),  b >= -a-sqrt(a^2-3), FALSE)
  eqn4 <- ifelse(abs(a+b) <= eps, b >= sqrt(3), FALSE)
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

  if(is.vector(moments))
    moments <- t(moments)
  
  # [a, b] = support of the standardized distribution
  sd <- sqrt(moments[, 2]-moments[, 1]^2) # strandard derivation
  a <- (lower-moments[, 1])/sd
  b <- (upper-moments[, 1])/sd
  g1 <- (moments[, 3]-3*sd^2*moments[, 1]-moments[, 1]^3)/sd^3
  g2 <- (-3*moments[, 1]^4+6*moments[, 2]*moments[, 1]^2-
           4*moments[, 3]*moments[, 1]+moments[, 4])/sd^4-3

  Q <- 4*g1*(a+b)+(3-a^2)*(3-b^2)
  
  # Existance of a distribution
  eqn1 <- ifelse(1+a*b <= 0, TRUE, FALSE)
  eqn2 <- ifelse(g2 >= g1^2-2, TRUE, FALSE)
  eqn3 <- ifelse(g2 <= g1^2-2-(a*b*(g1-a+1/a)*(g1-b+1/b))/(1+a*b), TRUE, FALSE)
  eqn4 <- ifelse(g1 <= b-1/b, TRUE, FALSE)
  eqn5 <- ifelse(g1 >= a-1/a, TRUE, FALSE)
  eqn6 <- (eqn1 & eqn2 & eqn3 & eqn4 & eqn5)
  
  ### Upper bound (Teuscher, Guiard)
  form1 <- 4*g1*(a^2+b^2+a*b-3)
  form2 <- ifelse(a+b > eps & Q>=0, 
                (form1-2*Q*(3+a*b+sqrt(Q))/(a+b))/(5*(a+b))-6/5, 0)
  form3 <- ifelse(a+b > eps & Q<0,
                (form1+Q*(4*g1+(a+b)*(a^2-3))/(3+2*a*b+a^2))/(5*(a+b))-6/5, 0)
  form4 <- ifelse(a+b < -eps & Q>=0,
                (form1-2*Q*(3+a*b-sqrt(Q))/(a+b))/(5*(a+b))-6/5, 0)
  form5 <- ifelse(a+b < -eps & Q<0,
                (form1+Q*(4*g1+(a+b)*(b^2-3))/(3+2*a*b+b^2))/(5*(a+b))-6/5, 0)
  form6 <- ifelse(abs(a+b) < eps, 
                1/5*((12*g1^2)/(3-a^2)+3*a^2-15), 0)
  formula <- form2 + form3 + form4 + form5 + form6
  eqn7 <- ifelse(g2 <= eps + formula, TRUE, FALSE)
 
  ### Lower bound (Johnsson, Rogers)
  r <- 2048*g1^2*(g1^2+sqrt(g1^2*(g1^2+4)))
  s <- sqrt(r^(1/3) - 256*g1^2*r^(-1/3))
  C <- 3 + s/2 - sqrt(128*g1^2/s - s^2)/2
  
  formula <- NULL
  if(g1 > eps) formula <- 6/5*(g1*(-sqrt(C)+1/sqrt(C))-1)
  if(g1 < -eps) formula <- 6/5*(g1*(sqrt(C)-1/sqrt(C))-1)
  if(abs(g1) < eps) formula <- -6/5
  
  eqn8 <- ifelse(g2 + eps >= formula, TRUE, FALSE)

  ### Combine everything
  
  results <- dplyr::case_when(
    isTRUE(eqn8 & eqn7) ~ "4-b-unimodal",
    isTRUE(!eqn6) ~ "not existant",
    TRUE ~ rep("not unimodal", length(eqn1))
  )
  
  return(results)
}