#======== todo =================================================================

#' Calculate raw moments
#'
#' This method calculates the raw moment of a given order for every time over
#' all simulations.
#' @param x object of class \code{\link{multSimCsv}}.
#' @param order number that specifies the order of the moments
#' @return data.frame with calculated moments
#' @seealso \code{\link[pdmpsim]{moments}} to calculate the moments of
#'   \code{\link{multSim}} objects
#' @importFrom LaF colmoment
#' @importFrom pdmpsim moments
#' @importFrom simecol fromtoby
#' @export
moments.multSimCsv <- function(x, order){
  if(!exists('colmoment', where = asNamespace('LaF'), mode = 'function')) 
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
  moments <- data.frame(order = order, time = times, moments)
  return(moments)
}