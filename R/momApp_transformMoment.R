##### transformMoment ####

#' A help function for \code{\link{momApp}}
#'
#' Let X be a random variable and Y the corresponding centered variable, i.e. \eqn{Y = X - \mu}. 
#' Let p be the \code{order} of the desired moment. 
#' If \code{type = 'raw'}, this function returns
#' \deqn{E(X^p) = \sum_{k_1=0}^{p_1}...\sum_{k_n=0}^{p_n} {p \choose k} \mu^{p-k} E(Y^k).}
#' The values of \eqn{E(Y^k)} are replaced by
#' \itemize{
#'   \item \code{centralMoment[[i]]} if \code{k} = \code{centralMomentOrders[i, ]},
#'   \item \code{\link{transformMoment}(type = 'central', ...)} if \code{k} is a row in \code{rawMomentOrders},
#'   \item \code{\link{symbolicMoments}(closure, k)} if \code{k} is neither a row in \code{centralMomentOrders} nor in \code{rawMomentOrders}.
#' }
#' If \code{type = 'central'}, this function returns
#' \deqn{E(Y^p) = \sum_{k_1=0}^{p_1}...\sum_{k_n=0}^{p_n} (-1)^{p-k} {p \choose k} \mu^{p-k} E(X^k).}
#' The values of \eqn{E(X^k)} are replaced by
#' \itemize{
#'   \item \code{rawMoment[[i]]} if \code{k} = \code{rawMomentOrders[i, ]},
#'   \item \code{\link{transformMoment}(type = 'raw', ...)} if \code{k} is not a row in \code{rawMomentOrders}.
#' }
#' 
#' @value This function returns an object of class 'momentList' with elements
#' \itemize{
#'   \item \code{centralMomentOrders}:
#'   \item centralMoments
#'   \item rawMomentOrders
#'   \item rawMoments
#'   \item closure
#' }
#' @param type string, either 'central' or 'raw'
#' @param order numeric vector giving the order of the desired moment
#' @param centralMomentOrders matrix. Every row gives the order of a central moment that is already known.
#' @param centralMoments list. The i-th entry is the central Moment of order \code{centralMomentOrders[i, ]}.
#' @param rawMomentOrders matrix. Every row gives the order of a raw moment that is already known.
#' @param rawMoments list. The i-th entry is the raw Moment of order \code{rawMomentOrders[i, ]}.
#' @param closure string giving the closure method to use if a central moment is unknown.
transformMoment <- function(order, type, momentList, closure = "zero"){
  
  momentList <- validate_momentList(momentList)
  
  p <- order
  n <- length(p)
  otherType <- setdiff(c("central", "raw"), type) # if type = 'raw' then otherType = 'central' and vice versa
  k_indexes <- as.matrix(expand.grid(lapply(1:n, function(i) 0:p[i])))

  readMoments <- function(momentList, type){ # a help function to assign some values
    
    if(type == "raw"){
      assign("typeOrders", momentList$rawMomentOrders, envir = parent.frame())
      assign("typeMoments", momentList$rawMoments, envir = parent.frame())
      assign("otherOrders", momentList$centralMomentOrders, envir = parent.frame())
      assign("otherMoments", momentList$centralMoments, envir = parent.frame())
    }
    if(type == "central"){
      assign("typeOrders", momentList$centralMomentOrders, envir = parent.frame())
      assign("typeMoments", momentList$centralMoments, envir = parent.frame())
      assign("otherOrders", momentList$rawMomentOrders, envir = parent.frame())
      assign("otherMoments", momentList$rawMoments, envir = parent.frame())
    }
    
    k_in_typeOrders <- apply(k_indexes, 1, function(i){
      if(is.null(typeOrders)) NA
      else prodlim::row.match(i, typeOrders)
    })
    assign("k_in_typeOrders", k_in_typeOrders, envir = parent.frame())
    
    k_in_otherOrders <- apply(k_indexes, 1, function(i){
      if(is.null(otherOrders)) NA
      else prodlim::row.match(i, otherOrders)
    })
    assign("k_in_otherOrders", k_in_otherOrders, envir = parent.frame())
  }
  readMoments(momentList, type)
  
  # before calculation: check if calculation is needed
  if(!is.na(prodlim::row.match(p, typeOrders))) # if there is already an entry in momentList
    return(momentList)
  if(type == "central" & is.na(prodlim::row.match(p, otherOrders))){
    moment <- symbolicMoments(distribution = closure, 
                              missingOrders = t(p),
                              mean = rep(0, n),
                              cov = cov.momentList(momentList))[[1]]
    
    momentList$centralMomentOrders <- rbind(momentList$centralMomentOrders, p)
    momentList$centralMoments <- append(momentList$centralMoments, moment)
    return(momentList)
  }
  
  # mu = vector with expected values
  muRowIndex <- lapply(1:n, function(i){
    muRow <- rep(0, n)
    muRow[i] <- 1
    prodlim::row.match(muRow, momentList$rawMomentOrders)
  })
  mu <- lapply(muRowIndex, function(i) momentList$rawMoments[[i]])
  if(sum(sapply(mu, is.null)) > 0){
    stop("momentList$rawMoments should contain all expected values.")
  }
  
  #### calculation #####
    
  sum <- list()
  for(k_index in 1:nrow(k_indexes)){
    
    k <- k_indexes[k_index, ]
    
    # if k is neither a row in typeMomentOrders nor in otherMomentOrders
    if(is.na(k_in_typeOrders[k_index]) & is.na(k_in_otherOrders[k_index])){
      momentCentral <- symbolicMoments(distribution = closure, 
                                       missingOrders = t(k),
                                       mean = mean(momentList),
                                       cov = cov.momentList(momentList))[[1]]
      momentList$centralMomentOrders <- rbind(momentList$centralMomentOrders, k)
      momentList$centralMoments <- append(momentList$centralMoments, momentCentral)
      readMoments(momentList, type)
    }
    # if k is a row in typeMomentOrders but not in otherMomentOrders
    if(!is.na(k_in_typeOrders[k_index]) & is.na(k_in_otherOrders[k_index])){
      momentList <- transformMoment(order = k, 
                                    type = otherType, 
                                    momentList = momentList, 
                                    closure = closure)
      readMoments(momentList, type)
    }
    # if k is a row in otherMomentOrders
    if(!is.na(k_in_otherOrders[k_index])){
      momentK <- otherMoments[[k_in_otherOrders[k_index]]]
    }
    
    # compute summand
    pChooseK <- as.numeric(Reduce("*", lapply(1:n, function(i) choose(p[i], k[i]))))
    muPower <- lapply(1:n, function(i) bquote(.(mu[[i]])^.(as.numeric(p[1:n]-k)[i])))
    muPower <- Reduce(function(a,b) bquote(.(a)*.(b)), muPower)
    if(type == "raw"){
      summand <- bquote(.(pChooseK)*.(muPower)*(.(momentK)))
    }
    if(type == "central"){
      sign <- (-1)^(sum(p)-sum(k))
      summand <- bquote(.(sign*pChooseK)*.(muPower)*(.(momentK)))
    }
    sum[[k_index]] <- summand
  }
  sum <- Reduce(function(a,b) bquote(.(a)+.(b)), sum)
  
  if(type == 'raw'){
    momentList$rawMomentOrders <- rbind(typeOrders, p)
    momentList$rawMoments <- append(typeMoments, sum)
  }
  if(type == 'central'){
    momentList$centralMomentOrders <- rbind(typeOrders, p)
    momentList$centralMoments <- append(typeMoments, sum)
  }
  return(momentList)
}  