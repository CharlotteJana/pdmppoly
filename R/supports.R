#======== todo =================================================================
#s1 PROM durch meine Promotion ersetzen
#s1 plotSupport ist eher sinnlos

#' Return the support for gene regulation models
#' 
#' In PROM, the support of the following gene regualtion models has been 
#' calculated: \code{\link[geneK]{Model K}}, \code{\link[geneK2]{Model K2}}, 
#' \code{\link[geneF]{Model F}}, \code{\link[geneKF]{Model KF}}, 
#' \code{\link[geneBF]{Model BF}} and \code{\link[geneT]{Model T}}.
#' Method \code{getSupport} returns a data.frame with the lower and upper
#' bounds of the support for every time value, method \code{plotSupport}
#' adds these bounds as lines into a plot given as \code{ggplot} object.
#' @param model Object of class \code{pdmpModel} or \code{polyPdmp}.
#' This has to be one of the previous mentioned models.
#' @examples
#' data(genePolyT)
#' getSupport(genePolyT) %>% head
#' 
#' data(genePdmpBF)
#' ms <-  multSim(genePdmpBF, seeds = 1:4)
#' plotSupport(ggplot = plotSeeds(ms, 1:4), data = getSupport(genePdmpBF))
#' @name support
#' @importFrom stringr str_detect
#' @aliases supports getSupport plotSupport getsupport plotsupport
getSupport <- function(model){
  t <- fromtoby(model@times)
  data <- data.frame(times = t)
  data <- with(as.list(c(model@parms, model@init)),{
    
    if(str_detect(descr(model), "Model K:")){
      data[, "lower"] <- f*exp(-b*t)
      data[, "upper"] <- f*exp(-b*t) + a/b*(1-exp(-b*t))
      return(data)
    }
    if(str_detect(descr(model), "Model K2:")){
      data[, "lower1"] <- f1*exp(-b1*t)
      data[, "upper1"] <- f1*exp(-b1*t) + a1/b1*(1-exp(-b1*t))
      data[, "lower2"] <- f2*exp(-b2*t) + (exp(-b2*t) - exp(-b1*t))*(f1*a2)/(b2-b1)
      data[, "upper2"] <- data$lower2 + (a1*a2)/(b2-b1)*(b1*(1-exp(-b2*t))-b2*(1-exp(-b1*t)))
      return(data)
    }
    if(str_detect(descr(model), "Model F:")){
      data[, "lower"] <- f*exp(-b*t)
      data[, "upper"] <- f*exp(-b*t) + a/b*(1-exp(-b*t))
      return(data)
    }
    if(str_detect(descr(model), "Model KF:")){
      data[, "lower"] <- f*exp(-b*t)
      data[, "upper"] <- f*exp(-b*t) + a/b*(1-exp(-b*t))
      return(data)
    }
    if(str_detect(descr(model), "Model BF:")){
      data[, "lower"] <- f*exp(-b*t) + a0/b*(1-exp(-b*t))
      data[, "upper"] <- f*exp(-b*t) + a1/b*(1-exp(-b*t))
      return(data)
    }
    if(str_detect(descr(model), "Model T:")){
      data[, "lower1"] <- fA*exp(-bA*t)
      data[, "lower2"] <- fB*exp(-bB*t)
      data[, "upper1"] <- fA*exp(-bA*t) + aA/bA*(1-exp(-bA*t))
      data[, "upper2"] <- fB*exp(-bB*t) + aB/bB*(1-exp(-bB*t))
      return(data)
    }
    stop("The support for this model is not included in the package. ", 
          "You have to calculate it yourself.")
  })
  return(data)
}

#' @rdname support
#' @param ggplot Object of class ggplot. The bounds will be plotted onto this 
#' object. If NULL, a new plot object will be created.
#' @param support data.frame that has the same structure as the data.frame
#' returned by \code{getSupport}.
#' @note Method \code{plotSupport} currently only works for models with one
#' continous variable.
plotSupport <- function(ggplot = NULL, support){
  if(is.null(ggplot))
    ggplot <- ggplot2::ggplot(data = NULL)
  
  plot <- ggplot + 
    ggplot2::geom_line(data = data, ggplot2::aes(x = times, y = lower), color = "blue") +
    ggplot2::geom_line(data = data, ggplot2::aes(x = times, y = upper), color = "blue")
  
  return(plot)
}