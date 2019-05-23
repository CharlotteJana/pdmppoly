#======== todo =================================================================
#t2: example in documentation
#t1: useCsv = multSimCsv implementieren
#t1: auskommentierten code mit manipulate bearbeiten
#t2: ms mit übergeben, so dass nur modality berechnet wird ODER
#    modality berechnungen auslagern
#t3: model und polyModel vergleichen
#t1: modality tests funktionieren nicht! 

#' Analyse a polynomial PDMP
#' 
#' Perform all methods that are available for a polynomial PDMP:
#' \itemize{
#' \item Simulate the PDMP multiple times
#' \item Plot the simulated data with all available plot methods
#' \item Calculate statistics such as mean, sd, median, etc
#' \item Perform moment approximation with all closure methods that are available
#' and plot the result
#' \item Test if the distribution is unimodal and plot the result
#' }
#' All plots and calculated data will be saved in different files.
#' It is possible to only perform a part of the methods, i.e. pass already simulated
#' data and only perform the plot methods or simulate without calculating statistics
#' or approximated moments.
#' 
#' @param polyModel an object of class \code{\link{polyPdmpModel}}
#' @param model the corresponding object of class \code{\link{pdmpModel}}
#' @param seeds number of seeds to be simulated
#' @param dir string giving the directory where files will be stored.
#' @param filenameprefix string. The name of each file saved by \code{analyseModel}
#' will start with this string.
#' @param momentorder integer giving the maximal order of moments that shall 
#' be calculated with \code{\link{momApp}}. Defaults to 10.
#' @param plotorder vector giving all the orders of moments that shall be plotted
#' with \code{\link{plotCompareMomApp}}. Defaults to 1:4.
#' @param useCsv boolean variable. Should \code{analyseModel} use 
#' \code{\link{multSim}} (if \code{FALSE}) or \code{\link{multSimCsv}} 
#' (if \code{TRUE}) for simulation? Defaults to \code{FALSE}, because
#' function \code{multSimCsv} is only necessary for huge simulations that 
#' exceed the working storage.
#' @param statistics character vector. Each entry should be the name of a function
#' that can be applied over simulated data. It will be used by the functions
#' \code{\link[pdmpsim]{plotStats}} and \code{\link[pdmpsim]{summarise_at}}.
#' @param plot boolean variable. Should \code{analyseModel} generate plots?
#' @param sim boolean variable. Should \code{analyseModel} do the simulation
#' or use already simulated data? In the latter case, it will use files
#' stored at the same places where \code{analyseModel} would store the simulations.
#' Furthermore, no value for parameter \code{seeds} is necessary.
#' @param momApp boolean variable. Should \code{analyseModel} calculate moments
#' with \code{\link{momApp}}?
#' @param modality boolean variable. Should \code{analyseModel} test if the 
#' distribution is unimodal with \code{\link{is.unimodal}}?
#' @param lower integer. Lower bound of the support of the distribution.
#' This variable is only needed if \code{modality = TRUE}.
#' @param upper integer. Upper bound of the support of the distribution.
#' This variable is only needed if \code{modality = TRUE}.
#' @importFrom pdmpsim format multSim discStates getMultSimData descr
#' @importFrom ggplot2 labs
#' @importFrom simecol "times<-" "init<-" "init" "parms"
#' @importFrom grDevices dev.off dev.print png
#' @importFrom momcalc is.unimodal
#' @export
analyseModel <- function(polyModel, model = polyModel, seeds = NULL, 
                         dir = file.path(getwd(), "simulations"), 
                         filenameprefix = descr(polyModel),
                         momentorder = 4, plotorder = 1:4, 
                         plot = TRUE, modality = TRUE, sim = TRUE, momApp = TRUE,
                         lower = NULL, upper = NULL, useCsv = FALSE,
                         statistics = c("min", "max", "mean", "median", "sd")){
  
  #### variables ####
  initNames <- names(init(model))
  parmsNames <- names(parms(model))
  discVars <- names(discStates(model))
  contVars <- setdiff(initNames, discVars)
  closureMethods <- c("zero", "zero", "normal", "lognormal", "gamma")
  closureCentral <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
  fname <- file.path(dir, filenameprefix)
  
  # to avoid the R CMD Check NOTE 'no visible binding for global variable ...'
  variable <- method <- time <- E <- . <- NULL
  
  message("\n", pdmpsim::format(model, short = F, collapse = "\n",
                                slots = c("descr", "parms", "init", "times")))
    
  #### simulation ####
  if(useCsv){
    
  }
  else{
    if(!sim){ # load existing simulated data
      
      ms <- readRDS(file = paste0(fname, ".rda"))
      message("Get MultSimData")
      msData <- readRDS(file = paste0(fname, "__multSimData.rda"))
      moments <- readRDS(file = paste0(fname, "__moments.rda"))
      
      if(!all.equal(model, ms$model))
        stop("Simulation stored in ", paste0(fname, ".rda"),
             " was done with another model or other model parameters.")
      if(is.null(seeds))
        seeds <- ms$seeds
      if(!identical(seeds, ms$seeds))
        stop("Simulation stored in ", paste0(fname, ".rda"),
             " was done with different seeds.")
    }
    else{
      ### simulation
      ms <- multSim(model, seeds)
      saveRDS(ms, file = paste0(fname, ".rda"))
      
      message("Get MultSimData")
      try({
        msData <- getMultSimData(ms)
        saveRDS(msData, file = paste0(fname, "__multSimData.rda"))
      })
    }
  }
  
  ### statistics 
  try({
    stats <- pdmpsim::summarise_at(
      msData,
      .vars = initNames,
      .funs = statistics)
    saveRDS(stats, paste0(fname,"__statistics.rda"))
  })
    
  #### moments ####
  if(momApp){
    try({
      message("Approximate Moments")
      ma <- momApp(polyModel, maxOrder = max(momentorder),
                   closure = closureMethods, centralize = closureCentral)
      ma <- addSimulation(ma, ms)
      saveRDS(ma, file = paste0(fname, "__moments.rda"))
    })
  }
  
  #### modality ####
    
  if(modality){
    message("Modality tests")
    
    # evaluate or set values for 'lower' and 'upper'
    values <- subset(msData, variable == name, select = value)
    if(is.null(lower)) 
      lower <- min(values)
    else
      lower <- with(as.list(parms(model)), eval(lower))
    if(is.null(upper)) 
      upper <- max(values)
    else
      upper <- with(as.list(parms(model)), eval(upper))
    
    modalities <- modalityTest(ma, lower, upper)
    saveRDS(modalities, file = paste0(fname, "__modality.rda"))
  }
    
  ##### plots #####
    
  if(plot & exists("msData")){
    
    # violin plot
    try({
      message("Plots: violin plot, ", appendLF = FALSE)
      pdmpsim::plotTimes(msData,
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
    })
    
    # boxplot with seednumbers
    try({
      message("boxplot, ", appendLF = FALSE)
      pdmpsim::plotTimes(msData,
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
    })
    
    # statistics
    try({
      message("statistics, ", appendLF = FALSE)
      pdmpsim::plotStats(msData,
                         vars = initNames[!(initNames %in% discVars)],
                         funs = statistics)
      ggplot2::ggsave(filename = paste0(fname,"__statistics.png"), 
                      dpi = 300, width = 20.4, height = 11, units = "cm")
    })
    
    # histogram for last simulated time value
    try({
      message("histogram, ", appendLF = FALSE)
      h <- hist(msData, t = times(model)["to"],
                main = descr(model),
                sub = pdmpsim::format(model, short = F, slots = "parms"))
      ggplot2::ggsave(filename = paste0(fname,"__histogram.png"), plot = h, 
                      dpi = 300, width = 20.4, height = 11, units = "cm")
      
      dev.off()
    })
    
    # densities for different time values
    try({
      message("densities, ", appendLF = FALSE)
      times <- unique(msData$time)
      times <- times[seq(1, length(times), length.out = 6)]
      times <- times[2:6]
      density(msData, t = times,
              main = descr(model),
              sub = pdmpsim::format(model, short = F, slots = "parms"))
      dev.print(png, filename = paste0(fname, "__densities.png"),
                width = 20.4, height = 11, units = "cm", res = 140)
      dev.off()
    })
    
    # moments
    try({
      message("moments, ", appendLF = FALSE)
      mplot <- plot(ma)

      ggplot2::ggsave(filename = paste0(fname,"__moments.png"), plot = mplot, 
                      dpi = 300, width = 20.4, height = length(model@init)*5.5, 
                      units = "cm")
    })
    
    # modality
      try({
        message("modality, ")
        timedist <- times(model)["by"]
        
        ggplot(data = modalities, aes(x = time, y = method, color = modality)) + 
          ggplot2::geom_segment(aes(xend = time + timedist, yend = method), 
                                size = 10, lineend = "butt") +
          ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=5))) +
          ggplot2::labs(
            title = model@descr,
            subtitle = format(model, slots = c("parms"), short = FALSE),
            caption = paste("Number of Simulations:", length(seeds))) +
          ggplot2::facet_wrap(~ variable, ncol = 1)
        
        ggplot2::ggsave(filename = paste0(fname,"__modality.png"), 
                        dpi = 300, width = 20.4, height = length(model@init)*5.5, 
                        units = "cm")
      })
    
    # plot (histogram over all simulations and times)
    if(!useCsv){
      try({
        message("overview ", appendLF = FALSE)
        plot(ms, discPlot = "line") # if you set discPlot = "smooth", you should suppress messages of ggsave
        ggplot2::ggsave(paste0(fname,"__plot.png"), dpi = 300, 
                        width = 20.4, height = 11, units = "cm")
      })
    }
  }

    
  #### end  ####
  message("All files are stored in ", dir)
}
