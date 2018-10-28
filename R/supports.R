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
      data[, "lower"] <- ξ*exp(-β*t)
      data[, "upper"] <- ξ*exp(-β*t) + α/β*(1-exp(-β*t))
      return(data)
    }
    if(str_detect(descr(model), "Model K2:")){
      data[, "lower1"] <- ξ1*exp(-β1*t)
      data[, "upper1"] <- ξ1*exp(-β1*t) + α1/β1*(1-exp(-β1*t))
      data[, "lower2"] <- ξ2*exp(-β2*t) + (exp(-β2*t) - exp(-β1*t))*(ξ1*α2)/(β2-β1)
      data[, "upper2"] <- data$lower2 + (α1*α2)/(β2-β1)*(β1*(1-exp(-β2*t))-β2*(1-exp(-β1*t)))
      return(data)
    }
    if(str_detect(descr(model), "Model F:")){
      data[, "lower"] <- ξ*exp(-β*t)
      data[, "upper"] <- ξ*exp(-β*t) + α/β*(1-exp(-β*t))
      return(data)
    }
    if(str_detect(descr(model), "Model KF:")){
      data[, "lower"] <- ξ*exp(-β*t)
      data[, "upper"] <- ξ*exp(-β*t) + α/β*(1-exp(-β*t))
      return(data)
    }
    if(str_detect(descr(model), "Model BF:")){
      data[, "lower"] <- ξ*exp(-β*t) + α0/β*(1-exp(-β*t))
      data[, "upper"] <- ξ*exp(-β*t) + α1/β*(1-exp(-β*t))
      return(data)
    }
    if(str_detect(descr(model), "Model T:")){
      data[, "lower1"] <- ξA*exp(-βA*t)
      data[, "lower2"] <- ξB*exp(-βB*t)
      data[, "upper1"] <- ξA*exp(-βA*t) + αA/βA*(1-exp(-βA*t))
      data[, "upper2"] <- ξB*exp(-βB*t) + αB/βB*(1-exp(-βB*t))
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
    ggplot2::geom_line(data = data, ggplot2::aes(x = times, y = lower), color = "red") +
    ggplot2::geom_line(data = data, ggplot2::aes(x = times, y = upper), color = "red")
  
  print(plot)
  return(invisible(plot))
}