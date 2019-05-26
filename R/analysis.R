#======== todo =================================================================
#t2: model und polyModel aus data auslesen
#t1: funktioniert das ganze? dir, subdirs überprüfen!
#t1: documentation
#t1: Funktion umbenennen

#' Analyse different polynomial PDMPs at once
#' 
#' Read a csv file which contains informations about the PDMPs
#' and their simulation parameters. Analyse all PDMPs with
#' \code{\link{analyseModel}} and save the resulting files in a 
#' given directory.
#' 
#' @param data a data.frame, see Details below
#' @param model an object of class \code{\link{pdmpModel}}
#' @param polyModel the corresponding object of class \code{\link{polyPdmpModel}}
#' @param dir string giving the directory where files will be stored.
#' @param subDirs boolean variable indicating if subdirectories for every
#' simulation shall be created. Defaults to FALSE.
#' @param ... additional parameters to function \code{\link{analyseModel}}.
#' @export
analysis <- function(data, model, polyModel, prefix, dir, subdirs, ...){
  
  for(i in seq_len(nrow(data))){
    
    #### filenames ####
    if(subDirs)
      dir.create(file.path(dir, data[i, "prefix"]), 
                 showWarnings = FALSE, recursive = TRUE)
    
    fname <- ifelse(subDirs, 
                    file.path(dir, data[i, "prefix"], data[i, "prefix"]),
                    file.path(dir, data[i, "prefix"]))
    
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
    
    analyseModel(model, polymodel, filenameprefix = fname, ...)
  }
}
