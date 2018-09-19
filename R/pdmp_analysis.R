#======== todo =================================================================
#t1: documentation
#t1: useCsv = multSimCsv implementieren
#t2: model und polyModel aus data auslesen
# install_github("CharlotteJana/pdmpsim@charlotte")

#' Analysis of models used in thesis ...
#' 
#' @param data a data.frame, see Details below
#' @param model an object of class \code{\link{pdmpModel}}
#' @param seeds number of seeds to be simulated
#' @param useCsv boolean variable indicating if \code{\link{multSim}} (FALSE)
#' or \code{\link{multSimCsv}} should be used for simulation. Defaults to FALSE.
#' @param dir string giving the directory where files will be stored.
#' @param subDirs boolean variable indicating if subdirectories for every
#' simulation shall be created. Defaults to FALSE.
#' @param momentorders vector giving all the orders of moments that shall 
#' be calculated. Defaults to 1:10.
#' @param plotorder vector giving all the orders of moments that shall be plotted.
#' Defaults to 1:4.
#' @importFrom pdmpsim format
#' @importFrom ggplot2 labs
#' @export
analysis <- function(data, model, polyModel, seeds = 1:50, useCsv = FALSE, 
                     dir = file.path(getwd(), "simulations"), subDirs = FALSE, 
                     momentorders = 1:10, plotorder = 1:4){
  
  #workingdir <- getwd()
  #setwd(dir)
  #on.exit(setwd(workingdir))
  
  initNames <- names(init(model))
  parmsNames <- names(parms(model))
  discVars <- names(discStates(model))
  
  for(i in seq_len(nrow(data))){
    ### set new values for init, parms, times ###
    suppressWarnings({
      for(name in initNames){
        try({
          if(!is.na(data[i, name])){
            init(model)[name] <- data[i, name]
            init(polyModel)[name] <- data[i, name]
          }
        }, silent = TRUE)
      }
      for(name in parmsNames){
        try({
          if(!is.na(data[i, name])){
            parms(model)[name] <- data[i, name]
            parms(polyModel)[name] <- data[i, name]
          }
        }, silent = TRUE)
      }
      for(name in c("from", "to", "by")){
        try({
          if(!is.na(data[i, name])){
            times(model)[name] <- data[i, name]
            times(polyModel)[name] <- data[i, name]
          }
        }, silent = TRUE)
      }
    })
  
    # filenames
    if(subDirs)
      dir.create(file.path(dir, data[i, "prefix"]), showWarnings = FALSE, recursive = TRUE)
    fname <- ifelse(subDirs, 
                    file.path(dir, data[i, "prefix"], data[i, "prefix"]),
                    file.path(dir, data[i, "prefix"]))
            #paste0(data[i,"prefix"], "_", 
            #pdmpsim::format(model, short = T, slots = "parms"))
    
    message("\n", pdmpsim::format(model, short = F, collapse = "\n",
                                  slots = c("descr", "parms", "init", "times")))
    ### simulate ###
    
    if(useCsv){
      
    }
    else{
      ms <- multSim(model, seeds)
      saveRDS(ms, file = paste0(fname, ".rda"))
      msData <- getMultSimData(ms)
      
      ### moments
      moments <- list()
      msim <- NULL
      for(m in momentorders){
        msim <- dplyr::bind_rows(msim, moments(ms, m))
      }
      moments[["Simulation"]] <- msim
      #saveRDS(moments, file = paste0(fname, "__moments.rda"))
      
      ### plot (histogram over all simulations and times)
      message("Plots: overview, ", appendLF = FALSE)
      suppressWarnings(plot(ms))
      ggplot2::ggsave(paste0(fname,"__plot.png"), dpi = 300, 
                      width = 20.4, height = 11, units = "cm")
     
    }
    
    ##### Plots #####
    
    # violin plot
    message("violin plot, ", appendLF = FALSE)
    plotTimes(msData,
              vars = initNames,
              plottype = "violin") +
    ggplot2::labs(title = descr(model),
                  subtitle = paste0("Number of simulations: ",
                                    length(unique(msData$seed)), "\n",
                                    pdmpsim::format(model, short = F,
                                                    collapse = "\n",
                                                    slots = "parms")))
    ggplot2::ggsave(filename = paste0(fname,"__violins.png"), dpi = 300,
                    width = 20.4, height = 11, units = "cm")


    # boxplot with seednumbers
    message("boxplot, ", appendLF = FALSE)
    plotTimes(msData,
              vars = initNames[!(initNames %in% discVars)],
              nolo = 3,
              plottype = "boxplot") +
      ggplot2::labs(title = descr(model),
                    subtitle = paste0("Number of simulations: ",
                                      length(unique(msData$seed)), "\n",
                                      pdmpsim::format(model, short = F,
                                                      collapse = "\n",
                                                      slots = "parms")))
    ggplot2::ggsave(filename = paste0(fname,"__boxplot.png"), dpi = 300,
                    width = 20.4, height = 11, units = "cm")

    # statistics (min, max, mean, median)
    #png(paste0(dir,"/", fname,"__statistics.png"), width = 1200, res = 140)
    message("statistics, ", appendLF = FALSE)
    plotStats(msData,
              vars = initNames[!(initNames %in% discVars)],
              funs = c("min", "mean", "median", "max"))
    ggplot2::ggsave(filename = paste0(fname,"__statistics.png"), 
                    dpi = 300, width = 20.4, height = 11, units = "cm")

    statistics <- summarise_at(msData,
                               .vars = initNames,
                               .funs = c("min", "max", "mean", "median", "sd"))
    saveRDS(statistics, paste0(fname,"__statistics.rda"))

    # histogram for last simulated time value
    message("histogram, ", appendLF = FALSE)
    h <- hist(msData, t = times(model)["to"],
              main = descr(model),
              sub = pdmpsim::format(model, short = F, slots = "parms"))
    ggplot2::ggsave(filename = paste0(fname,"__histogram.png"), plot = h, 
                    dpi = 300, width = 20.4, height = 11, units = "cm")

    # densities for different time values
    message("densities")
    times <- unique(msData$time)
    times <- times[seq(1, length(times), length.out = 6)]
    times <- times[2:6]
    dev.off()
    density(msData, t = times,
            main = descr(model),
            sub = pdmpsim::format(model, short = F, slots = "parms"))
    dev.print(png, filename = paste0(fname, "__densities.png"),
              width = 20.4, height = 11, units = "cm", res = 140)
    dev.off()
  
    #### moments ####
    
    message("Approximate Moments")
    for(s in c("reduceDegree", "setZero")){
      mcalc <- momApp(polyModel, l = max(momentorders), closure = s)$moments
      moments[[s]] <- mcalc
    }
    summary(moments)
    moments <- dplyr::bind_rows(moments, .id = "method")
    saveRDS(moments, file = paste0(fname, "__moments.rda"))
    
    message("Plot Moments")
    plotdata <- reshape2::melt(moments, 1:3, stringsAsFactors = TRUE)
    plotdata$method <- factor(plotdata$method, levels = c("Simulation", 
                                                          "reduceDegree", 
                                                          "setZero"))
    plotdata$variable <- as.character(plotdata$variable)
    plotdata <- subset(plotdata, order <= plotorder)

    plot <- ggplot2::ggplot(data = plotdata, ggplot2::aes(x = time, y = value))+ 
      ggplot2::geom_line(size = 1,
                         ggplot2::aes(color = method, linetype = method)) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(), 
                     axis.title.x = ggplot2::element_blank()) + 
      ggplot2::labs(
           title = model@descr,
           subtitle = format(model, slots = c("parms"), short = FALSE),
           caption = paste("Number of Simulations:", length(seeds))) +
      ggplot2::facet_wrap(variable ~ order, 
            scales = "free_y", nrow = length(model@init),
            labeller = ggplot2::label_bquote(cols = E(.(variable)^.(order))))
      ggplot2::ggsave(filename = paste0(fname,"__moments.png"), plot = plot, 
                      dpi = 300, width = 20.4, height = length(model@init)*5.5, 
                      units = "cm")

      message("All files are stored in ",ifelse(subDirs, file.path(dir, data[i, "prefix"]), dir))
  }
}