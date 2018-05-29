#======== todo =================================================================
#t1 LAF umdefinieren (mit Git!)
#t3 documentation nur für multSimCsv oder für beide?

#' Calculate raw moments
#'
#' This method calculates the raw moment of a given order for every time over
#' all simulations.
#' @param x object of class \code{\link{multSim}} or \code{\link{multSimCsv}}
#' @param order number that specifies the order of the moments
#' @return data.frame with calculated moments
#' @name moments
#' @export
moments <- function(x, order){
  UseMethod("moments", x)
}

#'@rdname moments
#'@importFrom LaF colmoment
#'@export
moments.multSimCsv <- function(x, order){
  if(!identical(find('colmoment'),"package:LaF")) 
    stop("Method 'colmoment' is not defined in package 'LaF'.")
  moments <- NULL
  times <- fromtoby(x$model@times)
  for(j in 1:length(x$lafList)){
    momentrow <- colmoment(x$lafList[[j]], 
                           columns = 2:(length(times)+1), 
                           order = order)
    moments <- cbind(moments, momentrow)
  }
  colnames(moments) <- sapply(1:length(x$csvList), function(i) 
    regmatches(x$csvList, regexec("_([^_]+).csv", x$csvList))[[i]][2])
  rownames(moments) <- NULL
  moments <- data.frame(time = times, moments)
  return(moments)
}