#======== todo =================================================================
#t2: example in documentation
#t3: simulierte daten direkt übergeben können
#t1: lower, upper mehrdimensional testen
#v1: Titel der modality plots
#t3: statistics nicht jedes mal berechnen
#t2: vergleich der modelle am anfang: mit kurzer zeit simulieren

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
#' @param model the corresponding object of class \code{\link{pdmpModel}}. It is
#'   used for simulation because the simulation of an object of class
#'   \code{pdmpModel} is faster than the simulation of an object of class
#'   \code{polyPdmpModel}.
#' @param seeds number of seeds to be simulated
#' @param dir string giving the directory where all files will be stored.
#' @param filenameprefix string. The name of each file saved by \code{analyseModel}
#' will start with this string.
#' @param momentorder numerical vector. Each entry specifies the maximal order
#'   of moments that shall be calculated with \code{\link{momApp}}. The moment
#'   approximation will be performed for each entry seperately. 
#'   Defaults to c(4, 10).
#' @param plotorder vector giving all the orders of moments that shall be plotted. 
#' Defaults to 1:4.
#  @param useCsv boolean variable. Should \code{analyseModel} use 
#  \code{\link{multSim}} (if \code{FALSE}) or \code{\link{multSimCsv}} 
#  (if \code{TRUE}) for simulation? Defaults to \code{FALSE}, because
#  function \code{multSimCsv} is only necessary for huge simulations that 
#  exceed the working storage.
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
#'@param lower numeric vector or matrix or data.frame specifying the lower
#'  bounds of the compact distribution that determines the law of the PDMP. #'
#'  It is an argument to function \code{\link{modalityTest}} and only needed if
#'  \code{modality = TRUE}. If \code{upper} . If \code{lower} is a vector, the
#'  i-th entry should give the lower bound of the i-th continous variable given
#'  in \code{init(model)}, independent of a time value. If \code{lower} is a
#'  matrix or data.frame, it should have a column named \code{time} containing
#'  all time values specified in \code{times(model)}. The other column names
#'  should be identical to the continous variables of the PDMP and contain the
#'  lower bounds for the corresponding variable and time value.
#' @param upper numeric vector or matrix or data.frame specifying the upper
#'   bounds of the compact distribution that determines the law of the PDMP. It
#'   is an argument to function \code{\link{modalityTest}} and only needed if
#'   \code{modality = TRUE}. If \code{upper} is a vector, the i-th entry should
#'   give the upper bound of the i-th continous variable given in
#'   \code{init(model)}, independent of a time value. If \code{upper} is a
#'   matrix or data.frame, it should have a column named \code{time} containing
#'   all time values specified in \code{times(model)}. The other column names
#'   should be identical to the continous variables of the PDMP and contain the
#'   upper bounds for the corresponding variable and time value.
#' @importFrom pdmpsim format multSim discStates getMultSimData descr
#' @importFrom ggplot2 labs aes ggplot
#' @importFrom simecol "times<-" "init<-" "init" "parms"
#' @importFrom grDevices dev.off dev.print png
#' @export
analyseModel <- function(polyModel, model = polyModel, seeds = NULL, 
                         dir = file.path(getwd(), "simulations"), 
                         filenameprefix = descr(polyModel),
                         momentorder = c(4,10), plotorder = 1:4, 
                         plot = TRUE, modality = TRUE, sim = TRUE, momApp = TRUE,
                         lower = NULL, upper = NULL, 
                         statistics = c("min", "max", "mean", "median", "sd"),
                         title = descr(model)){
  
  if(!identical(sim(polyModel, outSlot = FALSE, seed = 20),
                sim(model, outSlot = FALSE, seed = 20))){
    stop("Simulation of 'polyModel' and 'model' differ, 
          they do not represent the same PDMP.")
  }
  
  #### variables ####
  initNames <- names(init(model))
  parmsNames <- names(parms(model))
  discVars <- names(discStates(model))
  contVars <- setdiff(initNames, discVars)
  closureMethods <- c("zero", "zero", "normal", "lognormal", "gamma")
  closureCentral <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
  
  fname <- file.path(dir, filenameprefix)
  if(!dir.exists(dir)) dir.create(dir)
  
  # to avoid the R CMD Check NOTE 'no visible binding for global variable ...'
  variable <- method <- time <- E <- . <- NULL
  
  message("\n", pdmpsim::format(model, short = F, collapse = "\n",
                                slots = c("descr", "parms", "init", "times")))
    
  #### simulation ####
  if(!sim){ # load existing simulated data
    
    ms <- readRDS(file = paste0(fname, "__simulations.rda"))
    message("Get MultSimData")
    msData <- readRDS(file = paste0(fname, "__multSimData.rda"))
    
    if(!all.equal(model, ms$model))
      stop("Simulation stored in ", paste0(fname, "__simulatons.rda"),
           " was done with another model or other model parameters.")
    if(is.null(seeds))
      seeds <- ms$seeds
    if(!identical(seeds, ms$seeds))
      stop("Simulation stored in ", paste0(fname, "__simulations.rda"),
           " was done with different seeds.")
  }
  else{
    ### simulation
    ms <- multSim(model, seeds)
    saveRDS(ms, file = paste0(fname, "__simulations.rda"))
    
    message("Get MultSimData")
    try({
      msData <- getMultSimData(ms)
      saveRDS(msData, file = paste0(fname, "__multSimData.rda"))
    })
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
  ma <- list()
  for(i in seq_along(momentorder)){
    if(momApp){
      try({
        message("Approximate Moments of order <= ", momentorder[i])
        ma[[i]] <- momApp(polyModel, maxOrder = momentorder[i],
                     closure = closureMethods, centralize = closureCentral)
        ma[[i]] <- addSimulation(ma[[i]], ms)
        saveRDS(ma[[i]], file = paste0(fname, "__moments_order<=", momentorder[i], ".rda"))
      })
    }
    else{
      mafile <- paste0(fname, "__moments_order<=", momentorder[i], ".rda")
      if(file.exists(mafile))
        try({
          ma[[i]] <- readRDS(file = mafile)
        })
      else
        warning("file ", mafile, " doesn't exist.")
    }
  }
  
  #### modality ####
    
  if(modality){
    message("Modality tests")
    
    # set values for 'lower' and 'upper' if no values are provided
    setLower <- is.null(lower)
    setUpper <- is.null(upper)
    if(setLower | setUpper){
      for(i in seq_along(contVars)){
        values <- subset(msData, variable == contVars[i], select = value)
        if(setLower) 
          lower[i] <- min(values)
        if(setUpper) 
          upper[i] <- max(values)
      }
    }
    
    modalities <- list()
    for(i in seq_along(momentorder)){
      modalities[[i]] <- modalityTest(ma[[i]], lower, upper, vars = contVars)
      saveRDS(modalities[[i]], file = paste0(fname, "__modality_order<=", momentorder[i], ".rda"))
    }
  }
  else{
    modalities <- list()
    for(i in seq_along(momentorder)){
      modfile <- paste0(fname, "__modality_order<=", momentorder[i], ".rda")
      if(file.exists(modfile)){
        try({
          modalities[[i]] <- readRDS(file = modfile)
        })
      }
    }
  }
    
  ##### plots #####
    
  if(plot & exists("msData")){
    
    # violin plot
    try({
      message("Plots: violin plot, ", appendLF = FALSE)
      pdmpsim::plotTimes(msData,
                         vars = initNames,
                         plottype = "violin") +
        ggplot2::labs(title = title,
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
        ggplot2::labs(title = title,
                      subtitle = paste0("Number of simulations: ",
                                        length(unique(msData$seed)), "\n",
                                        pdmpsim::format(model, short = FALSE,
                                                        collapse = "\n",
                                                        slots = "parms")))
      ggplot2::ggsave(filename = paste0(fname,"__boxplot.png"), dpi = 300,
                      width = 20.4, height = 11, units = "cm")
    })
    
    # first seeds
    try({
      message("single simulations, ", appendLF = FALSE)
      pdmpsim::plotSeeds(msData, seeds = 1:4) +
        ggplot2::labs(title = title,
                      subtitle = pdmpsim::format(model, short = FALSE,
                                                 collapse = "\n",
                                                 slots = "parms"))
      ggplot2::ggsave(filename = paste0(fname,"__singleSimulations.png"), 
                      dpi = 300, width = 20.4, height = 11, units = "cm")
    })
    
    # statistics
    try({
      message("statistics, ", appendLF = FALSE)
      for(var in contVars){
        pdmpsim::plotStats(msData,
                           vars = var,
                           funs = statistics) +
        ggplot2::labs(title = title)
        ggplot2::ggsave(filename = paste0(fname,"__statistics_", var, ".png"), 
                        dpi = 300, width = 20.4, height = 11, units = "cm")
      }
    })
    
    # histogram for last simulated time value
    try({
      message("histogram, ", appendLF = FALSE)
      h <- hist(msData, t = times(model)["to"],
                main = title,
                sub = pdmpsim::format(model, short = F, slots = "parms"))
      dev.print(png, filename = paste0(fname, "__histogram.png"),
                width = 20.4, height = 11, units = "cm", res = 140)
      dev.off()
    })
    
    # densities for different time values
    try({
      message("densities, ", appendLF = FALSE)
      times <- unique(msData$time)
      times <- times[seq(1, length(times), length.out = 6)]
      times <- times[2:6]
      density(msData, t = times,
              main = title,
              sub = pdmpsim::format(model, short = F, slots = "parms"))
      dev.print(png, filename = paste0(fname, "__densities.png"),
                width = 20.4, height = 11, units = "cm", res = 140)
      dev.off()
    })
    
    # moments
    try({
      message("moments, ", appendLF = FALSE)
      for(i in seq_along(momentorder)){
        plot(ma[[i]]) + ggplot2::labs(title = paste0(title, " (closure for orders > ", momentorder[i],")"))
  
        ggplot2::ggsave(filename = paste0(fname, "__moments_order<=", momentorder[i], ".png"), 
                        dpi = 300, width = 20.4, height = length(model@init)*5.5, 
                        units = "cm")
      }
    })
    
    # modality
      try({
        message("modality, ", appendLF = FALSE)
        timedist <- times(model)["by"]
        
        for(i in seq_along(momentorder)){
          
          plotModalities(ma[[i]], modalities = modalities[[i]]) +
            ggplot2::labs(title = paste0(title, " (closure for orders > ", momentorder[i],")"))
          
          ggplot2::ggsave(filename = paste0(fname, "__modality_order<=", momentorder[i], ".png"), 
                          dpi = 300, width = 20.4, height = length(model@init)*5.5, 
                          units = "cm")
        }
      })
    
    # plot (histogram over all simulations and times)
    try({
      message("plot ")
      plot(ms, discPlot = "line") + # if you set discPlot = "smooth", you should suppress messages of ggsave
        ggplot2::labs(title = title)
      ggplot2::ggsave(paste0(fname,"__plot.png"), dpi = 300, 
                      width = 20.4, height = 11, units = "cm")
    })
  }

    
  #### end  ####
  message("All files are stored in ", dir)
}
