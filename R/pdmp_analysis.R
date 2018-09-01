#======== todo =================================================================
#t1: change threshold for boxplot
#t1: documentation
#t1: useCsv = multSimCsv implementieren

#' Analysis of models used in thesis ...
#' 
#' @param data a data.frame, see Details below
#' @param model an object of class \code{\link{pdmpModel}}
#' @param seeds number of seeds to be simulated
#' @param useCsv boolean variable, indicating if \code{\link{multSim}} (FALSE)
#' or \code{\link{multSimCsv}} should be used for simulation. Defaults to FALSE.
#' @param dir string giving the directory where files will be stored.
#' Defaults to the working directory.
#' @importFrom pdmpsim format
#' @importFrom ggplot2 labs
#' @export
analysis <- function(data, model, seeds = 1:50, useCsv = F, dir = getwd()){
  
  workingdir <- getwd()
  setwd(dir)
  on.exit(setwd(workingdir))
  
  initNames <- names(init(model))
  parmsNames <- names(parms(model))
  discVars <- names(discStates(model))
  
  message("All files are stored in ", dir)
  
  for(i in seq_len(nrow(data))){
    
    ### set new values for init, parms, times ###
    suppressWarnings({
      for(name in initNames){
        try({
          if(!is.na(data[i, name]))
            init(model)[name] <- data[i, name]
        }, silent = TRUE)
      }
      for(name in parmsNames){
        try({
          if(!is.na(data[i, name]))
            parms(model)[name] <- data[i, name]
        }, silent = TRUE)
      }
      for(name in c("from", "to", "by")){
        try({
        if(!is.na(data[i, name]))
          times(model)[name] <- data[i, name]
        }, silent = TRUE)
      }
    })
    message("\n", pdmpsim::format(model, short = F, collapse = "\n",
                                  slots = c("descr", "parms", "init", "times")))
    
    ### simulate ###
    
    if(useCsv){
      
    }
    else{
      ms <- multSim(model, seeds)
      fname <- paste0(data[i,1], "_", 
               pdmpsim::format(model, short = T, slots = "parms"))
      saveRDS(ms, file = paste0(fname, ".rda"))
      
      msData <- getMultSimData(ms)
      
      # plot histogram over all simulations and times
      plot(ms)
      ggplot2::ggsave(paste0(fname,"__plot.png"), dpi = 300, 
                      width = 20.4, height = 11, units = "cm")
    }
    
    # violin plot
    plotTimes(msData, 
              vars = initNames[!(initNames %in% discVars)],
              plottype = "violin") + 
    ggplot2::labs(title = descr(model),
                  subtitle = paste0("Number of simulations: ", 
                                    length(unique(msData$seed)), "\n",
                                    pdmpsim::format(model, short = F, 
                                                    collapse = "\n",
                                                    slots = "parms")))
    ggplot2::ggsave(filename = paste0(fname,"__violins.png"), dpi = 300)
    
    
    # boxplot with thresholds
    plotTimes(msData, 
              vars = initNames,
              threshold = 10,
              plottype = "boxplot") + 
      ggplot2::labs(title = descr(model),
                    subtitle = paste0("Number of simulations: ", 
                                      length(unique(msData$seed)), "\n",
                                      pdmpsim::format(model, short = F, 
                                                      collapse = "\n",
                                                      slots = "parms")))
    ggplot2::ggsave(filename = paste0(fname,"__boxplot.png"), dpi = 300)
    
    
    # histogram for last simulated time value
    hist(msData, t = times(model)["to"],
         main = descr(model),
         sub = pdmpsim::format(model, short = F, slots = "parms"))
    dev.print(png, filename = paste0(dir,"/", fname, "__histogram.png"), width = 1200, res = 140)
    dev.off()
    
    # densities for different time values
    times <- unique(msData$time)
    times <- times[seq(1, length(times), length.out = 6)]
    times <- times[2:6]
    density(msData, t = times,
         main = descr(model),
         sub = pdmpsim::format(model, short = F, slots = "parms"))
    dev.print(png, filename = paste0(dir,"/", fname, "__densities.png"), width = 1200, res = 140)
    dev.off()
    
    # statistics (min, max, mean, median)
    png(paste0(dir,"/", fname,"__statistics.png"), width = 1200, res = 140)
    plotStats(msData, 
              vars = initNames[!(initNames %in% discVars)],
              funs = c("min", "mean", "median", "max"))
    dev.off()
    
    statistics <- summarise_at(msData, 
                               .vars = initNames,
                               .funs = c("min", "max", "mean", "median", "sd"))
    saveRDS(statistics, paste0(fname,"__statistics.rda"))
  }
}