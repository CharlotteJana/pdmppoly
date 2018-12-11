checkParms <- function(model, parms = NULL, init = NULL, 
                       seeds = 5, data = NULL){
  if(!is.null(init)){
    for(name in names(init)){
      model@init[name] <- init[name]
    }
  }
  if(!is.null(parms)){
    for(name in names(parms)){
      model@parms[name] <- parms[name]
    }
  }
  #plot(sim(model, seed = seeds))
  ms <- multSim(model, seeds)
  plot <- plotSeeds(ms)
  # if(!is.null(data)){
  #   plot <- plot + 
  #     ggplot2::geom_line(data = data, ggplot2::aes(x = times, y = lower), color = "red") +
  #     ggplot2::geom_line(data = data, ggplot2::aes(x = times, y = upper))
  # }
  if(!is.null(data)){
    plot <- plotSupport(plot, data)
  }
  return(plot)
}

# library(manipulate)
# library(ggplot2)
# manipulate(checkParms(genePdmpKF, init = c(ξ = ξ), seeds = seeds,
#              parms = c(β = β, α = α, κ10 = κ10, κ01 = κ01, μ01 = μ01, μ10 = μ10)),
#            ξ = slider(0, 30),
#            β = slider(0.2, 10, step = 0.2),
#            α = slider(0.2, 10, step = 0.2),
#            κ10 = slider(0.2, 10, step = 0.2),
#            κ01 = slider(0.2, 10, step = 0.2),
#            μ10 = slider(0.2, 10, step = 0.2),
#            μ01 = slider(0.2, 10, step = 0.2),
#            seeds = slider(1,20))



# manipulate({
#   data <- data.frame(times = fromtoby(genePdmpBF@times))
#   data[, "lower"] <- ξ*exp(-β*data$times) + α0/β*(1-exp(-β*data$times))
#   data[, "upper"] <- ξ*exp(-β*data$times) + α1/β*(1-exp(-β*data$times))
#   plot <- checkParms(genePdmpBF, init = c(ξ = ξ), seeds = seeds,
#                       parms = c(β = β, α0 = α0, α1 = α1, κ10 = κ10, κ01 = κ01),
#                      data = data)
#   },
#            ξ = slider(0, 30),
#            β = slider(0.2, 10, step = 0.2),
#            α0 = slider(0.2, 10, step = 0.2),
#            α1 = slider(0.2, 10, step = 0.2),
#            κ10 = slider(0.2, 10, step = 0.2),
#            κ01 = slider(0.2, 10, step = 0.2),
#            seeds = slider(1,20))