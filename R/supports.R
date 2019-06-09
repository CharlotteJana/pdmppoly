#======== todo =================================================================
#t1 References: meine Promotion
#t2 Tests für getSupport

#' Return the support for gene regulation models
#' 
#' This package contains multiple PDMPs that simulate gene regulation
#' mechanisms. They can be loaded with \code{data(model)} where \code{model}
#' is one of the following: 
#' \itemize{
#' \item \code{\link{genePdmpK}}: gene regulation with constant activation
#' \item \code{\link{genePdmpK2}}: gene regulation with constant activation, 
#' transcription and translation are modeled seperately
#' \item \code{\link{genePdmpF}}: gene regulation with positive feedback
#' \item \code{\link{genePdmpBF}}: gene regulation with positive feedback
#' and basal transcription rate
#' \item \code{\link{genePdmpKF}}: gene regulation with constant activation
#' and positive feedback
#' \item \code{\link{toggleSwitch}}: toggle switch for two genes
#' } 
#' The distribution of each of these models has a compact support (see References
#' for further details). Method \code{getSupport} returns a list with two data.frames 
#' containing the lower and upper bounds of the support. They are time-depend 
#' and are given for every variable of the model seperately.
#' @param model Object of class \code{pdmpModel} or \code{polyPdmp}. This has to
#'   be one of the previous mentioned models.
#' @return A list with two entries named 'lower' and 'upper'. Both are data.frames
#' containing the values of the lower / upper bound for each variable of the PDMP
#' and each time value given in \code{times(model)}.
#' @examples
#' data(genePolyT)
#' getSupport(genePolyT)
#' 
#' @name support
#' @importFrom stringr str_detect
#' @aliases supports getSupport plotSupport getsupport plotsupport
getSupport <- function(model){
  t <- fromtoby(model@times)
  data <- list()
  data$lower <- data.frame(time = t)
  data$upper <- data.frame(time = t)
  data <- with(as.list(c(model@parms, model@init)),{
    
    if(str_detect(descr(model), "Model K:")){
      data$lower[, "f"] <- f*exp(-b*t)
      data$upper[, "f"] <- f*exp(-b*t) + a/b*(1-exp(-b*t))
      data$lower[, "d"] <- 0
      data$upper[, "d"] <- 1
      return(data)
    }
    if(str_detect(descr(model), "Model K2:")){
      data$lower[, "f1"] <- f1*exp(-b1*t)
      data$upper[, "f1"] <- f1*exp(-b1*t) + a1/b1*(1-exp(-b1*t))
      data$lower[, "f2"] <- f2*exp(-b2*t) + (exp(-b2*t) - exp(-b1*t))*(f1*a2)/(b2-b1)
      data$upper[, "f2"] <- data$lower[, "f2"] + (a1*a2)/(b1*b2*(b2-b1))*(b2*(1-exp(-b1*t))-b1*(1-exp(-b2*t)))
      data$lower[, "d"] <- 0
      data$upper[, "d"] <- 1
      return(data)
    }
    if(str_detect(descr(model), "Model F:")){
      data$lower[, "f"] <- f*exp(-b*t)
      data$upper[, "f"] <- f*exp(-b*t) + a/b*(1-exp(-b*t))
      data$lower[, "d"] <- 0
      data$upper[, "d"] <- 1
      return(data)
    }
    if(str_detect(descr(model), "Model KF:")){
      data$lower[, "f"] <- f*exp(-b*t)
      data$upper[, "f"] <- f*exp(-b*t) + a/b*(1-exp(-b*t))
      data$lower[, "d"] <- 0
      data$upper[, "d"] <- 1
      return(data)
    }
    if(str_detect(descr(model), "Model BF:")){
      data$lower[, "f"] <- f*exp(-b*t) + a0/b*(1-exp(-b*t))
      data$upper[, "f"] <- f*exp(-b*t) + a1/b*(1-exp(-b*t))
      data$lower[, "d"] <- 0
      data$upper[, "d"] <- 1
      return(data)
    }
    if(str_detect(descr(model), "Model T:")){
      data$lower[, "fA"] <- fA*exp(-bA*t)
      data$lower[, "fB"] <- fB*exp(-bB*t)
      data$upper[, "fA"] <- fA*exp(-bA*t) + aA/bA*(1-exp(-bA*t))
      data$upper[, "fB"] <- fB*exp(-bB*t) + aB/bB*(1-exp(-bB*t))
      data$lower[, "d"] <- 1
      data$upper[, "d"] <- 4
      return(data)
    }
    stop("The support for this model is not included in the package. ", 
          "You have to calculate it yourself.")
  })
  return(data)
}