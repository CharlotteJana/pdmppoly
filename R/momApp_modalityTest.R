#t1 lower und upper müssen vektoren sein, deren länge der anzahl von variablen entspricht

#' @importFrom tidyr spread
#' @importFrom momcalc is.unimodal
modalityTest <- function(momApp, lower, upper, 
                         vars = names(init(momApp$model))){
  
  
  stopifnot(length(lower) == length(vars))
  stopifnot(length(upper) == length(vars))
  
  modalities <- data.frame()
  
  for(calcMethod in unique(momApp$moments$method)){

    for(i in seq_along(vars)){
      
      # select moments
      m <- subset(momApp$moments, method == calcMethod & order <= 4, 
                  select = c("time", "order", vars[i]))
      m <- tidyr::spread(m, order, vars[i])
      m <- m[order(m$time),]
      m2 <- momcalc::is.unimodal(
        lower = lower[i], upper = upper[i],
        moments = m[, -1]
      )
      modalityMethod <- data.frame(time = m$time,
                                   method = calcMethod,
                                   variable = vars[i],
                                   modality = factor(m2, levels = c( "4-b-unimodal",
                                                                     "not unimodal",
                                                                     "not existant",
                                                                     NA_character_)))
      modalities <- dplyr::bind_rows(modalities, modalityMethod)
    }
  }
  modalities$method <- as.factor(modalities$method)
  # modalities <- tidyr::spread(modalities, variable, value)
  
  return(modalities)
}

#' @importFrom ggplot2 ggplot geom_segment guides labs facet_wrap
plotModalities <- function(momApp, modalities = NULL, 
                           vars = names(init(momApp$model)), ...){
  
  if(is.null(modalities))
    modalities <- modalityTest(momApp, ...)
  
  modalities <- subset(modalities, variable %in% vars)
  
  timedist <- times(momApp$model)["by"]
  
  plot <- ggplot(data = modalities, aes(x = time, y = method, color = modality)) + 
          geom_segment(aes(xend = time + timedist, yend = method), size = 5, lineend = "butt") +
          guides(colour = ggplot2::guide_legend(override.aes = list(size=5))) +
          labs(
            title = descr(momApp$model),
            subtitle = format(momApp$model, slots = c("parms"), short = FALSE)) +
          facet_wrap(~ variable, ncol = 1)
  
  if(!is.null(momApp$seeds))
    plot <- plot + 
            labs(caption = paste("number of simulations:", length(momApp$seeds)))
  
  return(plot)
}