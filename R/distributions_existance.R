#======== todo =================================================================
#t1 references for inequalities
#t3 exists.distribution -> integral = 1 prüfen?
#t3 messages umformulieren?
#t2 warum wurden nicht alle ungl aus Simpson_Welch_amNeuesten.wxm implementiert?

#' Check for existance of a distribution
#' 
#' This method checks different inequalities to determine if there exists a one
#' dimensional distribution with given support \eqn{[A, B]} and raw moments
#' \eqn{m_1, m_2, m_3, ...}{m₁, m₂, m₃, ...}. Only inequalities that include
#' moments up to order 5 are implemented.
#' @param A numeric. The lower bound of the support of the distribution.
#' @param B numeric. The upper bound of the support of the distriubtion.
#' @param m numeric vector giving the non standardized moments \eqn{m_1, m_2,
#'   m_3, ...}{m₁, m₂, m₃, ...}, sorted by their degree. 
#'   This vector should have at least two entries.
#' @examples
#' exists.distribution(0, 1, actuar::mbeta(1:4, 0.5, 0.5)) # TRUE
#' exists.distribution(-1, 2, actuar::munif(1:4, min = -1, max = 2)) # TRUE
#' exists.distribution(-4, -5, 1:5) # FALSE
#' @return boolean variable indicating if a distribution exists (TRUE) or 
#' not (FALSE). This method does not return a suitable distribution.
exists.distribution <- function(A, B, m){
  
  if(length(m) < 2){
    stop("At least two moment values are required.")
  }
  
  if(m[[2]] < m[[1]]^2)
    return(FALSE)
  
  #----- definitions ------
  
  bool <- TRUE
  sd <- sqrt(m[[2]]-m[[1]]^2) # standard derivation
  a <- (A-m[[1]])/sd # lower bound of support of the standardized distribution
  b <- (B-m[[1]])/sd # upper bount of support of the standardized distribution

  stand.moments <- NULL
  for(j in seq_along(m)){  # standardize moments m
    e <- m[[1]]
    s <- sapply(1:j, function(k) choose(j,k)*m[[k]]*(-e)^(j-k))
    stand.moments[[j]] <- (sum(s) + (-e)^j)/ sd^j
  }
  m <- stand.moments
  
  #----- inequalities for support ------
  
  # Integral = 1 prüfen?
  
  if(1+a*b > 0){ # 1.1, Teuscher, Guiard: (3)
    message("There is no distribution with lower bound ", A, 
            " and upper bound ", B)
    bool <- FALSE
  }
  
  #----- inequalities for 3rd moment ------
  
  if(length(m) >= 3){ 
    if(m[[3]] < a - 1/a){ # 1.2, Teuscher, Guiard: (3)
      message("There is no distribution with lower bound ", A,
              " and the 3rd moment given as ", m[[3]])
      bool <- FALSE
    }
    if(m[[3]] > b - 1/b){ # 1.3, Teuscher, Guiard: (3)
      message("There is no distribution with upper bound ", B,
              " and the 3rd moment given as ", m[[3]])
      bool <- FALSE
    }
  }
  
  #------- inequalities for 4rth moment ------
  
  if(length(m) >= 4){ # 1.4, Teuscher, Guiard: (2)
    if(m[4]-m[3]^2-1 > -(a^2-a*m[3]-1)*(b^2-b*m[3]-1) / (1+a*b)){
      message("There is no distribution with the 3rd moment being ", m[[3]],
              " the 4rth moment being ", m[[4]], " and the standardized ",
              "support being [", a, ",", b, "]")
      bool = FALSE
    }
    if(m[4]-m[3]^2-1 < 0){ # 2.6, Teuscher, Guiard: (1)
      message("There is no distribution with the 3rd moment being ", m[[3]],
              " and the 4rth moment being ", m[[4]])
      bool = FALSE
    }
  }
  
  #------- inequalities for 5th moment ------
  
  if(length(m) >= 5){ # 3.1
    if(-m[5] + 
       (m[4]*(m[4]-(b+2*a)*(m[3]+a^2*b)-b^2-2*a^2))/(m[3]-(a^2+1)*b-2*a) +
       (((b+a)^2+2*a^2)*m[3]^2)/(m[3]-(a^2+1)*b-2*a) + 
       ((a^3*b*(2*b+a)-2*a*(b+a)^2)*m[3])/(m[3]-(a^2+1)*b-2*a) + 
       (a^2*((b+a)^2+(2-a^2)*b^2))/(m[3]-(a^2+1)*b-2*a) < 0){
      message("There is no distribution with the given moments and boundaries")
      bool = FALSE
    }
    
    # 3.2
    r <- sqrt(a^2*m[4]^2 + 
                ((6*a-2*a^3)*m[3]-6*a^2+4)*m[4] - 
                4*a*m[3]^3 + 
                (a^4+6*a^2-3)*m[3]^2 + 
                (2*a-6*a^3)*m[3] + 
                4*a^4-3*a^2)/(2*a^2) +
      (a^2*(3-m[4]) + a^3*m[3] - 3*a*m[3] - 2)/(2*a^3)
    if(-(r^(1/3)* 
         (-r*a^7*m[5] + (r*a^8-4*r*a^6)*m[4] - 6*r*a^6*m[3]^2 +
          (16*r*a^7-18*r*a^5)*m[3] - 6*r*a^8 + 18*r*a^6 - 9*r*a^4) + 
         r^(2/3)*
         (4*r*a^7*m[4] - 2*a^5*m[3]^3 + (6*a^6-6*a^4)*m[3]^2 + 
          (-4*r*a^8-6*a^7+12*r*a^6+12*a^5-6*a^3)*m[3] + 
          (r^2+2)*a^8 - 12*r*a^7 - 6*a^6 + 8*r*a^5 + 6*a^4 -2*a^2) +
         (4*r*a^6*m[3]-4*r*a^7+4*r*a^5)*m[4] + 
         a^4*m[3]^4 + (4*a^3-4*a^5)*m[3]^3 + 
         (-4*r*a^7+6*a^6+12*r*a^5-12*a^4+6*a^2)*m[3]^2 + 
         (4*r*a^8+(-2*r^2-4)*a^7-28*r*a^6+12*a^5+20*r*a^4-12*a^3+4*a)*m[3] +
         (2*r^2+1)*a^8+12*r*a^7+(-2*r^2-4)*a^6-20*r*a^5+6*a^4 +8*r*a^3-4*a^2+1
    )/(r^(4/3)*a^7) < 0){
      message("There is no distribution with the given moments and boundaries")
      bool = FALSE
    }
  }
  
  return(bool)
}