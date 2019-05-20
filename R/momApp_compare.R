
#' @inheritParams momApp obj maxOrder
#' @param closures character vector
#' @param centralize logical vector
compareMomApp <- function(obj, maxOrder = 4, ms = NULL,
                          closures = c("zero", "zero", "normal", "lognormal", "gamma"),
                          centralize = c(TRUE, FALSE, TRUE, FALSE, FALSE)){
  moments <- list()
  methodNames <- closures
  
  for(s in seq_along(closures)){
    
    if(closures[s] == "normal")
      methodNames[s] <- "normal (central)"
    if(closures[s] == "zero" & centralize[s])
      methodNames[s] <- "zero (central)"
    if(closures[s] == "zero" & !centralize[s])
      methodNames[s] <- "zero (raw)"
    
    try({
      mcalc <- momApp(obj = obj, 
                      maxOrder = maxOrder, 
                      closure = closures[s],
                      centralize = centralize[s])$moments
      mcalc <- cbind(method = methodNames[s],
                     mcalc)
      moments[[s]] <- mcalc
    })
  }
  if(!is.null(ms)){
    msim <- NULL
    for(m in seq_len(maxOrder)){
      msim <- dplyr::bind_rows(msim, moments(ms, m))
    }
    msim <- cbind(method = "simulation",
                  msim)
    moments[[length(closures) + 1]] <- msim
  }
  moments <- dplyr::bind_rows(moments)
  moments$method <- as.factor(moments$method)
  return(moments)
}

plotCompareMomApp <- function(model, moments, maxOrder = 4, simnumber = NA_character_){
  
  plotdata <- reshape2::melt(moments, 1:3, stringsAsFactors = TRUE)
  #plotdata$method <- factor(plotdata$method, levels = c("Simulation", 
  #                                                      approxMethods))
  plotdata$variable <- as.character(plotdata$variable)
  plotdata <- subset(plotdata, order <= maxOrder)
  
  mplot <- ggplot(data = plotdata, ggplot2::aes(x = time, y = value))+ 
    ggplot2::geom_line(size = 1,
                       ggplot2::aes(color = method, linetype = method)) +
    ggplot2::scale_linetype_manual(values = c(2, 3, 4, 1, 5, 6, 7)) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(), 
                   axis.title.x = ggplot2::element_blank()) +
    ggplot2::labs(
      title = model@descr,
      subtitle = format(model, slots = c("parms"), short = FALSE),
      caption = paste("Number of Simulations:", simnumber)) +
    ggplot2::facet_wrap(variable ~ order,
                        scales = "free_y", nrow = length(model@init),
                        labeller = ggplot2::label_bquote(cols = E(.(variable)^.(order))))

  return(mplot)
}