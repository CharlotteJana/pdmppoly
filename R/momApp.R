#======== todo =================================================================
#t1 tests überarbeiten
#t1 documentation von closure anpassen

#' Moment approximation for polynomial PDMPs
#' 
#' @param obj object of class \code{\link{polyPdmpModel}}.
#' @param maxOrder integer defining the highest order of moments that are considered.
#'   Higher degrees are droped and replaced by other values. The replacement
#'   method is specified in parameter \code{closure}.
#' @param closure string defining the method that does the moment closure, i. e.
#'   that changes the system of ODEs into a closed form that is solvable.
#'   Possible values are \code{zero} and ... .
#' @details 
#' The returned s3 class \code{momApp} contains 6 different elements:
#' \itemize{
#' \item \code{model}: the polyPdmpModel \code{obj}
#' \item \code{out}: a matrix giving the result of function \code{\link[deSolve]{ode}}.
#' Column names were added to make the result understandable.
#' \item \code{maxOrder}: integer defining the highest order of moments to be calculated
#' \item \code{moments}: a data.frame with the resulting moments, of the same structure
#'  as the result of function \code{\link[pdmpsim]{moments}} in package \pkg{pdmpsim}.
#' }
#' @examples 
#' data(genePolyK2)
#' a <- momApp(genePolyK2, maxOrder = 4)
#' plot(a)
#' print(a)
#' summary(a)
#' @name momApp
#' @aliases momentApproximation momentapproximation momapp
#' @include polypdmp_class.R polypdmp_accessors.R polypdmp_generator.R
#' @seealso \code{\link{momApp-methods}} for further analysing the result.
#' @return an object of class \code{momApp} (see details).
#' @note This method works only for PDMPs with one discrete variable.
#' @importFrom prodlim row.match
#' @importFrom spray index
#' @importFrom simecol fromtoby
#' @importFrom deSolve ode
#' @importFrom dplyr %>%
#' @importFrom momcalc momentList transformMoment symbolicMoments extractCov extractMean
#' @export
setGeneric("momApp",
           function(obj, maxOrder = 4, closure = "zero", 
                    centralize = TRUE, na.rm = TRUE)
             standardGeneric("momApp"))

#' @rdname momApp
#' @export
setMethod("momApp", signature(obj = "polyPdmpModel"), 
          function(obj, maxOrder = 4, closure = "zero", 
                   centralize = TRUE, na.rm = TRUE) {
  
  states <- obj@discStates[[1]]
  n <- length(obj@init) - length(obj@discStates) # continuous variables
  k <- length(states)
  names <- names(obj@init)  # names of all variables
  dname <- names(obj@discStates)
  cnames <- names[!names %in% dname]
  
  # to avoid the R CMD Check NOTE 'no visible binding for global variable ...'
  # time <- variable <- NULL
    
  ###### create all moment combinations that are needed ######
  
    # r = all moment combinations that are needed
    r <- data.frame(sapply(1:n, function(i) i = 0:maxOrder))
    r <- expand.grid(r)
    if(n > 1) r <- r[which(rowSums(r) >= 0 & rowSums(r[,1:n]) <= maxOrder), ]
    dimnames(r) <- list(1:nrow(r), cnames)
    
    # create k indicator variables that indicate the state of the discrete var
    indicatorNames <- sapply(1:k, function(i) paste0(dname, states[i]))
    indicatorMatrix <- matrix(rep(t(diag(k)), length.out = k*k*nrow(r)), 
                              ncol = k, byrow = TRUE)
   
    # lhs = all moment combinations times all indicator variables 
    #       (this will be the lefthandside of the ode system)
    lhs <- as.data.frame(r[rep(seq_len(nrow(r)), each=k),])
    lhs <- cbind(lhs, indicatorMatrix)
    dimnames(lhs) <- list(1:nrow(lhs), c(cnames, indicatorNames))
    lhs <- as.matrix(lhs)
    
  ###### create system of odes #####
    
    #calculate EVGenerator to every moment combination 
    odeList <- lapply(1:nrow(lhs), function(c)
      EVGenerator(obj, m = lhs[c,1:n], j = states[which(lhs[c, (n+1):(n+k)] == 1)]))
    
    matchingRows <- lapply(odeList, function(ode){
      sapply(1:nrow(index(ode)), 
             function(i) prodlim::row.match(index(ode)[i,], lhs))
    })
    
  ###### check for missing ode's and eventually perform moment closure ######
    
    # rowsToChange = index of ode's that contain moments which are not in lhs
    rowsToChange <- which(sapply(1:length(matchingRows), function(row) anyNA(matchingRows[[row]])))
    
    if(length(rowsToChange) == 0){
      lhsFull <- lhs
    }
    else{
      
      # lhsMissing = moments which are not in lhs (each row represents one moment, as in lhs)
      lhsMissing <- NULL
      for(i in rowsToChange){
        h <- which(is.na(matchingRows[[i]]))
        lhsMissing <- rbind(lhsMissing, index(odeList[[i]])[h, ])
      }
      
      colnames(lhsMissing) <- colnames(lhs)
      rownames(lhsMissing) <- rep("missing", nrow(lhsMissing))
      lhsMissing <- unique(lhsMissing)
      
      lhsFull <- rbind(lhs, lhsMissing)
      matchingRows <- lapply(odeList, function(ode){
        sapply(1:nrow(index(ode)), 
               function(i) prodlim::row.match(index(ode)[i,], lhsFull))
      })
    }
    
  ###### convert ode's into quoted formulas #######
    
    if(length(rowsToChange > 0)){ 
      
      missingMoments <- list()

      for(i in 1:nrow(lhsMissing)){
        missingRow <- lhsMissing[i, ]
        indicatorIndex <- which(missingRow[(n+1):(n+k)] == 1)
        indicatorLhs <- lhs[k*(1:(nrow(lhs)/k)-1) + indicatorIndex, 1:n, drop = FALSE]
        stateList <- lapply(rownames(indicatorLhs), function(name) bquote(state[.(as.numeric(name))]))
        mList <- momcalc::momentList(rawMomentOrders = indicatorLhs, 
                                     rawMoments = stateList, 
                                     warnings = FALSE)
        
        if(centralize){ # centralize ode's which contain missing moments and perform moment closure
          if(closure %in% c("lognormal", "gamma")){
            stop("Closure method '", closure, "' is only implemented for raw moments.")
          }

          mListExpanded <- momcalc::transformMoment(order = missingRow[1:n],
                                                    type = "raw",
                                                    closure = closure,
                                                    momentList = mList)
          
          missingMoments[[i]] <- mListExpanded$rawMoments[[prodlim::row.match(missingRow[1:n], mListExpanded$rawMomentOrders)]]
        }
        else{ # perform moment closure directly
          if(closure %in% c("normal")){
            stop("Closure method '", closure, "' is only implemented for centralized moments.")
          }
          if(closure == "zero"){
            cov <- NA
            mean <- NA
          }
          else{
            cov <- momcalc::extractCov(mList)
            mean <- momcalc::extractMean(mList)
          }
          missingMoments[[i]] <- momcalc::symbolicMoments(distribution = closure,
                                                 missingOrders = missingRow[1:n],
                                                 cov = cov, mean = mean)[[1]]
        }
      }
    }

    odeSystem <- rep(list(NA), nrow(lhs))
    for(j in seq_len(nrow(lhs))){
      list <- lapply(seq_along(matchingRows[[j]]), function(i){
        momentIndex <- matchingRows[[j]][[i]]
        momentIsMissing <- identical(rownames(lhsFull)[momentIndex], "missing")
        if(momentIsMissing){
          h <-  prodlim::row.match(lhsFull[momentIndex, ], lhsMissing)
          return(bquote(.(value(odeList[[j]])[i])*(.(missingMoments[[h]]))))
        }
        else{
          return(bquote(.(value(odeList[[j]])[i])*state[.(momentIndex)]))
        }
      })
      odeSystem[[j]] <- Reduce(function(a,b) bquote(.(a)+.(b)), list)
    }
    
  #### simulate the system of odes with deSolve ######
    
    discInd <- getIndex(obj@init[dname], states)
    state <- apply(lhs, 1, function(row) 
        1/length(discStates(obj)[[1]])*Reduce("*", obj@init[1:n]^row[1:n]) 
    ) #initial value: dirac messure with peak in 
      # obj@init["f"] and every discrete state
    times <- fromtoby(obj@times)
    func <- function(lhs, state, parms){
      list(sapply(odeSystem, function(x) eval(x)))
    }
   out <- ode(y = state, times = times, func = func, 
              parms = obj@parms, method = obj@solver)
   
  ### create meaningful column names for out
   
   contNames <- c(sapply(2:(nrow(r)), function(row)
     paste(cnames, r[row, 1:n], sep = "^", collapse = "*")
   ))
   contNames <- gsub("\\*?[^\\*\\(]+\\^0\\*?", "", contNames) # remove ^0
   contNames <- gsub("+\\^+1", "", contNames) # remove ^1
   discNames <- sapply(1:k, function(i) paste0("P(", dname, "=", states[i],")"))
   outNames <- NULL
   for(i in 1:k){
     outNames <- rbind(outNames, paste0("E(", contNames, "|", 
                                       dname, "=", states[i], ")"))
   }
   outNames <- c("time", discNames, outNames)
   colnames(out) <- outNames

  #### create data.frame with raw moments #####
    
   moments <- expand.grid(time = times, order = 1:maxOrder)
   moments[, names] <- NA
   
   # moments of discrete variables
   colnames <- paste0("P(", dname, "=", states, ")")
   for(j in 1:maxOrder){
    values <- rowSums(out[, colnames] %*% diag(states))
    moments[which(moments$order == j), dname] <- values
   }
   
   # moments of order 1
   for(i in 1:n){
     colnames <- paste0("E(", cnames[i], "|", dname, "=", states, ")")
     moments[which(moments$order == 1), cnames[i]] <- rowSums(out[, colnames])
   }
   
   # moments of order > 1
   for(i in 1:n){
     for(j in 2:maxOrder){
      colnames <- paste0("E(", cnames[i], "^", j,"|", dname, "=", states, ")")
      moments[which(moments$order == j), cnames[i]] <- rowSums(out[, colnames])
     }
   }
   if(na.rm == TRUE)
     moments <- moments[!is.na(rowSums(moments)), ]
    
  #### create class 'momApp' #####

  result <- structure(list(model = obj, 
                           moments = moments,
                           out = out,
                           maxOrder = maxOrder,
                           closure = closure,
                           centralize = centralize
                           ), class = "momApp")
    
   return(result)
})
   