#t1 Überall cov.momentList durch cov ersetzen

#' @name momentList
momentList <- function(rawMomentOrders = NULL,
                       rawMoments = list(),
                       centralMomentOrders = NULL,
                       centralMoments = list()){
  
  mList <- structure(list(rawMomentOrders = rawMomentOrders,
                          rawMoments = rawMoments,
                          centralMomentOrders = centralMomentOrders,
                          centralMoments = centralMoments,
                     class = "momentList"))
  return(mList)
}

#' @describeIn momentList
validate_momentList <- function(x){
  
  stopifnot(is.list(x$rawMoments))
  stopifnot(is.list(x$centralMoments))
  stopifnot(is.matrix(x$rawMomentOrders) | is.data.frame(x$rawMomentOrders))
  stopifnot(is.matrix(x$centralMomentOrders) | is.data.frame(x$centralMomentOrders))
  
  if(length(unique(x$rawMomentOrders)) != length(x$rawMomentOrders)){
    stop("Some rows in 'rawMomentOrders' appear several times. 
         They should only appear once.")
  }
  if(length(unique(x$centralMomentOrders)) != length(x$centralMomentOrders)){
    stop("Some rows in 'centralMomentOrders' appear several times. 
         They should only appear once.")
  }
  if(nrow(x$rawMomentOrders) != length(x$rawMoments)){
    stop("The number of elements in 'rawMoments' should be equal to
         the number of rows in 'rawMomentOrders'.")
  }
  if(nrow(x$centralMomentOrders) != length(x$centralMoments)){
    stop("The number of elements in 'centralMoments' should be equal to
         the number of rows in 'centralMomentOrders'.")
  }
  if(ncol(x$rawMomentOrders) != ncol(x$centralMomentOrders)){
    stop("The number of columns in 'rawMomentOrders' and 'centralMomentOrders'
         should be identical.")
  }
  
  # ------- check moments of order 1 -------
  
  n <- ncol(x$rawMomentOrders)
  for(i in seq_len(n)){
    unitVector <- rep(0, n)
    unitVector[n-i+1] <- 1
    rowIndex <- prodlim::row.match(unitVector, x$centralMomentOrders)
    if(!is.na(rowIndex)){
      if(x$centralMoments[[rowIndex]] != 0) 
        warning("Central moments of order 1 should be 0.")
    }
    else{
      x$centralMomentOrders <- rbind(unitVector, x$centralMomentOrders)
      x$centralMoments <- append(0, x$centralMoments)
    }
  }
  rownames(x$centralMomentOrders) <- NULL
  
  #------- check moments of order 0 -------
  
  indexRaw <- prodlim::row.match(rep(0, n), x$rawMomentOrders)
  indexCentr <- prodlim::row.match(rep(0, n), x$centralMomentOrders)
  
  if(!is.na(indexRaw)){
    if(x$rawMoments[[indexRaw]] != 1)
      warning("Moments of order 0 should have value 1.")
  }
  if(!is.na(indexCentr)){
    if(x$centralMoments[[indexCentr]] != 1)
      warning("Moments of order 0 should have value 1.")
  }
  
  if(is.na(indexRaw)){
    x$rawMomentOrders <- rbind(rep(0, n), x$rawMomentOrders)
    x$rawMoments <- append(1, x$rawMoments)
  }
  if(is.na(indexCentr)){
    x$centralMomentOrders <- rbind(rep(0, n), x$centralMomentOrders)
    x$centralMoments <- append(1, x$centralMoments)
  }
  
  return(x)
}

#' @describeIn momentList
cov.momentList <- function(x){
  n <- ncol(x$centralMomentOrders)
  cov <- list()
  
  for(i in 1:n){
    cov[[i]] <- list()
    for(j in 1:i){
      row <- rep(0, n)
      row[i] <- row[i] + 1
      row[j] <- row[j] + 1

      cov[[i]][[j]] <- x$centralMoments[[prodlim::row.match(row, x$centralMomentOrders)]]
      cov[[j]][[i]] <- cov[[i]][[j]]
    }
  }
  
  return(cov)  
}

#' @describeIn momentList
mean.momentList <- function(x){
  n <- ncol(x$rawMomentOrders)
  mean <- list()
  
  for(i in 1:n){
    row <- rep(0, n)
    row[i] <- 1
    mean[[i]] <- x$rawMoments[[prodlim::row.match(row, x$rawMomentOrders)]]
  }
  
  return(mean)
}
