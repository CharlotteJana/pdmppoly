#======== todo =================================================================
#t3 benötigte Pakete: distr <- wofür?

context("is.unimodal")

test_that("is.unimodal returns TRUE for simple unimodal distributions", {
  
  # uniform distribution:
  expect_true(is.unimodal(-10, 20, actuar::munif(1:4, min = -10, max = 20)))
  
  # beta distribution:
  expect_true(is.unimodal(0, 1, actuar::mbeta(1:4, 1, 2)))
})


test_that("is.unimodal returns FALSE for bimodal distributions", {
   
    # beta distribution:
    expect_false(is.unimodal(0, 1, actuar::mbeta(1:4, 0.5, 0.5)))
    
    # arcsine on [0,20]: 
    # see http://www.randomservices.org/random/special/Arcsine.html
    expect_false(is.unimodal(0, 20, c(20^1*1/2, 
                                      20^2*(1*3)/(2*4), 
                                      20^3*(1*3*5)/(2*4*6), 
                                      20^4*(1*3*5*7)/(2*4*6*8))))
})

################## Mixtures #####################

## 2 Gleichverteilungen
# is.unimodal(0, 5, mmixunif(2))                            # ✔ (unimodal)
# is.unimodal(0, 6, mmixunif(2, a = c(1, 4), b = c(2, 5)))  # ✘  (bimodal)
# is.unimodal(0, 6, mmixunif(2, a = c(1, 4), b = c(2, 5), weights = c(.3,.7))) # ✔
# is.unimodal(0, 10, mmixunif(2, a = c(1, 8), b = c(2, 9))) # ✔✘? (bimodal)


## 2 Lognormalverteilungen:
# is.unimodal(0, 20, mmixlnorm(2))        # Fehlermeldung, weil Support zu klein
# is.unimodal(0, 10000000, mmixlnorm(2))  # Fehlermeldung, weil Support zu klein
# is.unimodal(0, 1e15, mmixlnorm(2))      # ✘ (bimodal)

############### Genregulation ###################


# var <- "ξ"
# parms <- c(α = 1, γ = 0.36, κ01 = 0.36, κ10 = 0.3)
# #m <- EW[grep(paste("^", var,"[1-9]?$", sep = ""), names(EW), perl=TRUE)]
# 
# EWξ_model1 <- function(parms){
#   EW <- c()
#   EW["ξ"]  = with(as.list(parms), 
#                   (α*κ01)/(γ*(κ10+κ01)))
#   EW["ξ2"] = with(as.list(parms), 
#                   (α^2*κ01*(κ01+γ))/(γ^2*(κ10+κ01)*(κ10+κ01+γ)))
#   EW["ξ3"] = with(as.list(parms), 
#                   (α^3*κ01*(κ01+γ)*(κ01+2*γ)) / 
#                     (γ^3*(κ10+κ01)*(κ10+κ01+γ)*(κ10+κ01+2*γ)))
#   EW["ξ4"] = with(as.list(parms), 
#                   (α^4*κ01*(κ01+γ)*(κ01+2*γ)*(κ01+3*γ)) / 
#                     (γ^4*(κ10+κ01)*(κ10+κ01+γ)*(κ10+κ01+2*γ)*(κ10+κ01+3*γ)))
#   return(EW)
# }
# 
# is.unimodal(0, parms[["α"]]/parms[["γ"]], EWξ_model1(parms))