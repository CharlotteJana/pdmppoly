#======== todo =================================================================
#t1: documentation
#t1: useCsv = multSimCsv implementieren
#t1: auskommentierten code mit manipulate bearbeiten
#t2: model und polyModel aus data auslesen
#s1: Prom durch meine Promotion ersetzen
#t2: modality plot for every variable
#t1: modality tests funktionieren nicht! (siehe plots)
#t1: lower und upper müssen vektoren sein!

#' Analysis of models used in PROM
#' 
#' @param data a data.frame, see Details below
#' @param model an object of class \code{\link{pdmpModel}}
#' @param polyModel the corresponding object of class \code{\link{polyPdmpModel}}
#' @param seeds number of seeds to be simulated
#' @param useCsv boolean variable indicating if \code{\link{multSim}} (if FALSE)
#' or \code{\link{multSimCsv}} (if TRUE) should be used for simulation. 
#' @param dir string giving the directory where files will be stored.
#' @param subDirs boolean variable indicating if subdirectories for every
#' simulation shall be created. Defaults to FALSE.
#' @param momentorders vector giving all the orders of moments that shall 
#' be calculated. Defaults to 1:10.
#' @param plotorder vector giving all the orders of moments that shall be plotted.
#' Defaults to 1:4.
#' @param plot boolean variable. Should \code{analysis} generate plots?
#' @param sim boolean variable. Should \code{analysis} do the simulation
#' or use already simulated data? In the latter case, it will use files
#' stored at the same places where \code{analysis} would store the simulations.
#' @param lower integer. Lower bound of the compact support of the distribution.
#' @param upper integer. Upper bound of the compact support of the distribution.
#' @importFrom pdmpsim format multSim discStates getMultSimData descr
#' @importFrom ggplot2 labs
#' @importFrom simecol "times<-" "init<-" "init" "parms"
#' @importFrom grDevices dev.off dev.print png
#' @export
analysis <- function(data, model, polyModel, seeds = NULL, useCsv = FALSE, 
                     dir = file.path(getwd(), "simulations"), subDirs = FALSE, 
                     momentorders = 1:10, plotorder = 1:4, plot = TRUE, modality = TRUE,
                     sim = TRUE, lower = NULL, upper = NULL){

  #### variables ####
  initNames <- names(init(model))
  parmsNames <- names(parms(model))
  discVars <- names(discStates(model))
  contVars <- setdiff(initNames, discVars)
  approxMethods <- c("setZero", "reduceDegree")
  
  # to avoid the R CMD Check NOTE 'no visible binding for global variable ...'
  variable <- method <- time <- E <- . <- NULL
  
  for(i in seq_len(nrow(data))){
    
    #### filenames ####
    if(subDirs)
      dir.create(file.path(dir, data[i, "prefix"]), 
                 showWarnings = FALSE, recursive = TRUE)
    
    fname <- ifelse(subDirs, 
                    file.path(dir, data[i, "prefix"], data[i, "prefix"]),
                    file.path(dir, data[i, "prefix"]))
    
    #con <- file(paste0(fname, "_messages.txt"), open = "wt")
    #sink(con, type = "message")
    
    #### set new values for init, parms, times ####
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
      if(is.null(seeds)){
        try({
          if(!is.na(data[i, "maxSeed"])){
            seeds <- 1:data[i, "maxSeed"]
          }
        })
      }
    })
    
    message("\n", pdmpsim::format(model, short = F, collapse = "\n",
                                  slots = c("descr", "parms", "init", "times")))
    
    #### simulation ####
    if(useCsv){
      
    }
    else{
      if(!sim){
        ms <- readRDS(file = paste0(fname, ".rda"))
        #message("Get MultSimData")
        msData <- readRDS(file = paste0(fname, "__multSimData.rda"))
        moments <- readRDS(file = paste0(fname, "__moments.rda"))
        
        if(!identical(model, ms$model))
          stop("Simulation stored in ", paste0(fname, "__moments.rda"),
               " was done with other model parameters.")
        if(!identical(seeds, ms$seeds))
          stop("Simulation stored in ", paste0(fname, "__moments.rda"),
               " was done with different seeds.")

      }
      else{
        ### simulation
        ms <- multSim(model, seeds)
        saveRDS(ms, file = paste0(fname, ".rda"))
        
        message("Get MultSimData")
        msData <- getMultSimData(ms)
        saveRDS(msData, file = paste0(fname, "__multSimData.rda"))
        
        ### statistics
        try({
          statistics <- pdmpsim::summarise_at(
                                msData,
                                .vars = initNames,
                                .funs = c("min", "max", "mean", "median", "sd"))
          saveRDS(statistics, paste0(fname,"__statistics.rda"))
        })
        
        ### moments
        try({
          message("Approximate Moments")
          
          moments <- list()
          msim <- NULL
          for(m in momentorders){
            msim <- dplyr::bind_rows(msim, moments(ms, m))
          }
          moments[["Simulation"]] <- msim
          
          for(s in approxMethods){
            mcalc <- momApp(polyModel, max(momentorders), closure = s)$moments
            moments[[s]] <- mcalc
          }
          moments <- dplyr::bind_rows(moments, .id = "method")
          moments$method <- as.factor(moments$method)
          saveRDS(moments, file = paste0(fname, "__moments.rda"))
        })
        
      }
     
    }
    
    #### modality ####
    
    if(modality){
      message("Modality tests")
      modalities <- data.frame()
      
      for(calcMethod in c("Simulation", approxMethods)){
        modalityMethod <- data.frame(time = fromtoby(times(model)))
        
        for(name in contVars){
          
          if(is.null(lower) | is.null(upper)){ # set values if no support is given
            values <- subset(msData, variable == name, select = value)
            if(is.null(lower)) lower <- min(values)
            if(is.null(upper)) upper <- max(values)
          }
          
          # select moments
          m <- subset(moments, method == calcMethod & order <= 4, 
                      select = c("time", "order", name))
          m <- tidyr::spread(m, order, name)
          m <- m[order(m$time),]
          m2 <- apply(within(m, rm("time")), 1, 
                      function(row){
                        is.unimodal(
                          lower = with(as.list(parms(model)), eval(lower)), 
                          upper = with(as.list(parms(model)), eval(upper)), 
                          moments = row)
                      })
          colname <- paste("modality of", name)
          modalityMethod[, colname] <- factor(m2, levels = c( "4-b-unimodal",
                                                              "not unimodal",
                                                              "not existant",
                                                              NA_character_))
          modalityMethod[, "method"] <- calcMethod
        }
        modalities <- dplyr::bind_rows(modalities, modalityMethod)
      }
      modalities$method <- as.factor(modalities$method)
      modalities$method
      saveRDS(modalities, file = paste0(fname, "__modality.rda"))
    }
    
    ##### plots #####
    
    if(plot){
      
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
      
      # statistics (min, max, mean, median)
      try({
        message("statistics, ", appendLF = FALSE)
        pdmpsim::plotStats(msData,
                           vars = initNames[!(initNames %in% discVars)],
                           funs = c("min", "mean", "median", "max"))
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
        plotdata <- reshape2::melt(moments, 1:3, stringsAsFactors = TRUE)
        plotdata$method <- factor(plotdata$method, levels = c("Simulation", 
                                                              approxMethods))
        plotdata$variable <- as.character(plotdata$variable)
        plotdata <- subset(plotdata, order <= plotorder)
        
        mplot <- ggplot(data = plotdata, ggplot2::aes(x = time, y = value))+ 
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
        
        ggplot2::ggsave(filename = paste0(fname,"__moments.png"), plot = mplot, 
                        dpi = 300, width = 20.4, height = length(model@init)*5.5, 
                        units = "cm")
      })
      
      # modality (für mehrere variablen anpassen!)
      if(modality){
        try({
          message("modality, ")
          timedist <- times(model)["by"]
          
          ggplot(data = modalities, aes(x = time, y = method, color = `modality of f`)) + 
            ggplot2::geom_segment(aes(xend = time + timedist, yend = method), 
                                      size = 20, lineend = "butt") +
            ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=5))) +
            ggplot2::labs(
              title = model@descr,
              subtitle = format(model, slots = c("parms"), short = FALSE),
              caption = paste("Number of Simulations:", length(seeds)))
          
          ggplot2::ggsave(filename = paste0(fname,"__modality.png"), 
                          dpi = 300, width = 20.4, height = length(model@init)*5.5, 
                          units = "cm")
        })
      }
      
      # plot (histogram over all simulations and times)
      if(!useCsv){
        try({
          message("overview ", appendLF = FALSE)
          suppressMessages(suppressWarnings(plot(ms)))
          ggplot2::ggsave(paste0(fname,"__plot.png"), dpi = 300, 
                          width = 20.4, height = 11, units = "cm")
        })
      }
    }
    
    #### end  ####
    message("All files are stored in ",ifelse(subDirs, 
                                              file.path(dir, data[i, "prefix"]), 
                                              dir))
    
    #sink(type = "message")
  }
}
