#t1 Überall cov.momentList durch cov ersetzen

#' @name momentList
momentList <- function(rawMomentOrders = NULL,
                       rawMoments = list(),
                       centralMomentOrders = NULL,
                       centralMoments = list()){
  
  if(is.null(rawMomentOrders) & is.null(centralMomentOrders))
    stop("Please provide either values for 'rawMomentOrders' and 'rawMoments'
          or values for 'centralMomentOrders' and 'centralMoments'")
  
  if(is.null(rawMomentOrders)){
    rawMomentOrders <- t(rep(0, ncol(centralMomentOrders)))
    rawMoments <- list(1)
  }
  if(is.null(centralMomentOrders)){
    centralMomentOrders <- rbind(rep(0, ncol(rawMomentOrders)), diag(ncol(rawMomentOrders)))
    centralMoments <- append(1, as.list(rep(0, ncol(rawMomentOrders))))
  }
  
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
  if(nrow(x$rawMomentOrders) != 0 &
     nrow(x$centralMomentOrders) != 0 & 
     ncol(x$rawMomentOrders) != ncol(x$centralMomentOrders)){
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
  
  n <- ifelse(length(x$centralMoments) > 0, 
              ncol(x$centralMomentOrders), 
              ncol(x$rawMomentOrders))

  covOrder <- expand.grid(lapply(1:n, function(i) 0:2))
  covOrder <- covOrder[which(rowSums(covOrder) == 2),]
  covOrderMissing <- is.na(prodlim::row.match(covOrder, x$centralMomentOrders))
  
  for(i in seq_along(covOrderMissing)){
    if(covOrderMissing[i]){
      x <- transformMoment(order = as.numeric(covOrder[i, ]),
                           type = "central",
                           momentList = x)  
    }
  }
  
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
