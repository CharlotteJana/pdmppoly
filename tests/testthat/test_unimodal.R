######## Testbeispiele (bekannte Verteilungen) ########## 

# benötigte Pakete: actuar, distr
A <- -10
B <- 20

######## Unimodale Verteilungen ###########

#is.unimodal(A, B, munif(1:4, min = A, max = B)) # ✔
#is.unimodal(0, 1, mbeta(1:4, 1, 2))             # ✔

## Die folgenden Verteilungen haben eigentlich keinen kompakten Support -> Problem???????
#
is.unimodal(0, 20, actuar::mexp(1:4, rate = 0.3))                 # ✔
is.unimodal(-10, 20, actuar::mnorm(1:4, mean = 0, sd = 1))         # ✔
is.unimodal(-10, 20, actuar::mnorm(1:4, mean = 99, sd = 1))         # hier gibt es keine Verteilung ???
#is.unimodal(0, B, mlnorm(1:4, meanlog = 0, sdlog = 0.3)) # ✔
#is.unimodal(0, B, mlnorm(1:4, meanlog = 0, sdlog = 1.5)) # hier gibt es keine Verteilung ???

########## Bimodale Verteilungen ###########

#is.unimodal(0, 1, mbeta(1:4, 0.5, 0.5)) # ✔
#is.unimodal(0, B, c(B/2, (3*B^2)/8, (5*B^3)/16, (35*B^4)/128)) # ✔ Arcsine auf [0,B], siehe http://www.math.uah.edu/stat/special/Arcsine.html

################## Mixtures #####################

## 2 Gleichverteilungen
# is.unimodal(0, 5, mmixunif(2))                            # ✔ (unimodal)
# is.unimodal(0, 6, mmixunif(2, a = c(1, 4), b = c(2, 5)))  # ✘  (bimodal)
# is.unimodal(0, 6, mmixunif(2, a = c(1, 4), b = c(2, 5), weights = c(.3,.7))) # ✔
# is.unimodal(0, 10, mmixunif(2, a = c(1, 8), b = c(2, 9))) # ✔✘? (bimodal)

## 2 Normalverteilungen:
# is.unimodal(0, 20, mmixnorm(2, means = c(3, 13), weights = c(.3, .7)))  # ✔
# is.unimodal(0, 20, mmixnorm(2, means = c(3, 13)))                       # ✘
# is.unimodal(0, 20, mmixnorm(2, weights = c(.3, .7)))                    # ✘
# is.unimodal(0, 20, mmixnorm(2))                                         # ✘

## 2 Lognormalverteilungen:
# is.unimodal(0, 20, mmixlnorm(2))        # Fehlermeldung, weil Support zu klein
# is.unimodal(0, 10000000, mmixlnorm(2))  # Fehlermeldung, weil Support zu klein
# is.unimodal(0, 1e15, mmixlnorm(2))      # ✘ (bimodal)

############### Genregulation ###################

var <- "ξ"
parms <- c(α = 1, γ = 0.36, κ01 = 0.36, κ10 = 0.3)
#m <- EW[grep(paste("^", var,"[1-9]?$", sep = ""), names(EW), perl=TRUE)]

EWξ_model1 <- function(parms){
  EW <- c()
  EW["ξ"]  = with(as.list(parms), (α*κ01)/(γ*(κ10+κ01)))
  EW["ξ2"] = with(as.list(parms), (α^2*κ01*(κ01+γ))/(γ^2*(κ10+κ01)*(κ10+κ01+γ)))
  EW["ξ3"] = with(as.list(parms), (α^3*κ01*(κ01+γ)*(κ01+2*γ))/(γ^3*(κ10+κ01)*(κ10+κ01+γ)*(κ10+κ01+2*γ)))
  EW["ξ4"] = with(as.list(parms), (α^4*κ01*(κ01+γ)*(κ01+2*γ)*(κ01+3*γ))/(γ^4*(κ10+κ01)*(κ10+κ01+γ)*(κ10+κ01+2*γ)*(κ10+κ01+3*γ)))
  return(EW)
}

is.unimodal(0, parms[["α"]]/parms[["γ"]], EWξ_model1(parms))
