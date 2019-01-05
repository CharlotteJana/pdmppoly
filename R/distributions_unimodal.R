#======== todo =================================================================
#v1 bei unterschiedlichen ergebnissen wie reagieren?
#   bsp: 2-b-unimodal + not existant?

#' Test if moments come from a unimodal distribution with compact support
#'
#' Given a compact support \code{[A, B]} and moments \code{m1, m2, ...} of a 
#' one dimensional distribution, this method checks if the (unknown) 
#' distribution can be unimodal. \cr 
#' There exist several inequalities that test for nonunimodality (see References). 
#' Depending on the inequality, moments up to order 2 or 4 are
#' required. A distribution that satisfies all inequalities that contain only
#' moments up to order 2 is called \emph{2-b-unimodal}. A distribution that
#' satisfies all inequalities that contain only moments up to order 4 is called
#' \emph{4-b-unimodal}. The internal methods \code{is.2_b_unimodal} and
#' \code{is.4_b_unimodal} test these inequalities. Method \code{is.unimodal}
#' performs these checks, depending on the number of given distributions. A
#' 4-b-unimodal distribution can still be nonunimodal (see the examples below).
#' But a failure of these tests assures the distribution not to be unimodal.
#'
#' @param lower numeric or vector. The lower bound(s) A of the support \eqn{[A,
#'   B]} of the distribution.
#' @param upper numeric or vector. The upper bound(s) B of the support \eqn{[A,
#'   B]} of the distribution.
#' @param moments numeric vector giving the non standardized moments \eqn{m_1,
#'   m_2, m_3, ...}{m₁, m₂, m₃, ...}, sorted by their degree. This vector should
#'   have at least two entries. It is also possible to enter a matrix, where
#'   every row contains the moments \eqn{m_1, m_2, m_3, ...}{m₁, m₂, m₃, ...}. 
#' @param eps numeric value. Some inequalities are of the form \code{... > 0}.
#'   For numerical reasons it is better to test for \code{... > eps} where
#'   \code{eps} is a small number.
#' @return Character vector giving the results of the test. Possible values are
#' "not unimodal", "not existant", NA_character_ and "2-b-unimodal" or 
#' "4-b-unimodal".
#' @example inst/examples/unimodal.R
#' @references 
#' \insertRef{TeuscherGuiard1994}{pdmppoly}
#' 
#' \insertRef{JohnsonRogers1951}{pdmppoly}
#' 
#' \insertRef{SimpsonWelch1960}{pdmppoly}
#' @name is.unimodal
#' @importFrom Rdpack reprompt
#' @export
is.unimodal <- function(lower, upper, moments, eps = 1e-10){
  #[A,B] = Support der ZG
  #m = sortierter Vektor mit Momenten
  
  l <- c(length(lower), length(upper), nrow(rbind(moments)))
  l <- l[! l %in% 1] # remove all values that are 1
  if(length(unique(l)) > 1){
    stop("length of input vectors and matrix do not fit")
  }
  
  if(ncol(rbind(moments)) < 4){
    r2 <- is.2_b_unimodal(lower, upper, moments, eps)
  }
  
  if(ncol(rbind(moments)) >= 4){
    r1 <- is.2_b_unimodal(lower, upper, moments, eps)
    r2 <- is.4_b_unimodal(lower, upper, moments, eps)
    
    index <- which(r1 != "2-b-unimodal")
    if(any(r1[index] != r2[index]))
      warning("Results of is.2_b_unimodal and is.4_b_unimodal differ.")
  } 
  return(r2)
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
    is.na(eqn1 & eqn5) ~ NA_character_,
    !eqn1 ~ "not existant",
    eqn5 ~ "2-b-unimodal",
    TRUE ~ rep("not unimodal", length(eqn1))
  )
  
  return(results)
}

#' @rdname is.unimodal
#' @export
is.4_b_unimodal <- function(lower, upper, moments, eps = 1e-10){

  if(is.vector(moments))
    moments <- t(moments)
  
  # [a, b] = support of the standardized distribution
  sd <- sqrt(moments[, 2]-moments[, 1]^2) # standard derivation
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
  
  form1 <- ifelse(g1 > eps & !is.na(g1), 6/5*(g1*(-sqrt(C)+1/sqrt(C))-1), 0)
  form2 <- ifelse(g1 < -eps & !is.na(g1), 6/5*(g1*(sqrt(C)-1/sqrt(C))-1), 0)
  form3 <- ifelse(abs(g1) < eps & !is.na(g1), -6/5, 0)
  
  eqn8 <- ifelse(!is.na(g2) & g2 + eps >= form1 + form2 + form3, TRUE, FALSE)
  eqn8[is.na(g1) | is.na(g2)] <- NA

  ### Combine everything
  results <- dplyr::case_when(
    is.na(eqn8 & eqn7 & eqn6) ~ NA_character_,
    !eqn6 ~ "not existant",
    eqn8 & eqn7 ~ "4-b-unimodal",
    TRUE ~ rep("not unimodal", length(eqn1))
  )
  
  return(results)
}