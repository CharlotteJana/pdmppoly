#======== todo =================================================================
#v3 contRes und discRes umbenennen
#t2 ist contInd wirklich nötig?

#' Moment approximation for polynomial PDMPs
#' 
#' @param obj object of class \code{\link{polyPdmpModel}}.
#' @param maxOrder integer defining the highest order of moments that are considered.
#'   Higher degrees are droped and replaced by other values. The replacement
#'   method is specified in parameter \code{closure}.
#' @param closure string defining the method that does the moment closure, i. e.
#'   that changes the system of ODEs into a closed form that is solvable.
#'   Possible values are \code{setZero} (the default) and reduceDegree.
#' @details 
#' The returned s3 class \code{momApp} contains 6 different elements:
#' \itemize{
#' \item \code{model}: the polyPdmpModel \code{obj}
#' \item \code{discRes}: a matrix giving the calculated moments of the different
#' indicator variables that replace the discrete variable (see \code{\link{blowupSpray}}
#' for an explanation of the indicator variables)
#' \item \code{contRes}: a matrix giving the calculated moments of the continous variables
#' \item \code{maxOrder}: integer defining the highest order of moments to be calculated
#' \item \code{closure}: string giving the closure method. See ... for more details.
#' \item \code{contInd}: a data.frame with all moment indexes that are calculated
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
#' @seealso \code{\link{momentClosure}} for the internal method that performs
#'   the moment closure, \code{\link{momApp-methods}} for further analysing
#'   the result.
#' @return an object of class \code{momApp} (see details).
#' @note This method works only for PDMPs with one discrete variable.
#' @importFrom prodlim row.match
#' @importFrom spray index
#' @importFrom simecol fromtoby
#' @importFrom deSolve ode
#' @importFrom stats aggregate as.formula
#' @importFrom dplyr %>%
#' @importFrom momcalc momentList transformMoment symbolicMoments extractCov extractMean
#' @export
setGeneric("momApp",
           function(obj, maxOrder = 4, closure = "zero", centralize = TRUE)
             standardGeneric("momApp"))

#' @rdname momApp
#' @export
setMethod("momApp", signature(obj = "polyPdmpModel"), 
          function(obj, maxOrder = 4, closure = "zero", centralize = TRUE) {
  
  states <- obj@discStates[[1]]
  n <- length(obj@init) - length(obj@discStates) # continuous variables
  k <- length(states)
  names <- names(obj@init)  # names of all variables
  dnames <- names(obj@discStates)
  cnames <- names[!names %in% dnames]
  
  # to avoid the R CMD Check NOTE 'no visible binding for global variable ...'
  time <- variable <- NULL
    
  ###### create all moment combinations that are needed ######
  
    # r = all moment combinations that are needed
    r <- data.frame(sapply(1:n, function(i) i = 0:maxOrder))
    r <- expand.grid(r)
    if(n > 1) r <- r[which(rowSums(r) >= 0 & rowSums(r[,1:n]) <= maxOrder), ]
    dimnames(r) <- list(1:nrow(r), cnames)
    
    # create k indicator variables that indicate the state of the discrete var
    indicatorNames <- sapply(1:k, function(i) paste0(dnames, states[i]))
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
    print(missingMoments)
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
    
    discInd <- getIndex(obj@init[dnames], states)
    state <- apply(lhs, 1, function(row) 
        1/length(discStates(obj)[[1]])*Reduce("*", obj@init[1:n]^row[1:n]) 
    ) #initial value: dirac messure with peak in obj@init["f"] and ...
    times <- fromtoby(obj@times)
    func <- function(lhs, state, parms){
      list(sapply(odeSystem, function(x) eval(x)))
    }
    out <- ode(y = state, times = times, func = func, 
               parms = obj@parms, method = obj@solver)

  #### create class 'momApp' from the result 'out' #####
    
    # discRes = only indicator variables
    discRes <- out[, c(1:(k+1))] # Achtung: Zeitspalte kommt dazu
    colnames(discRes)[-1] <- colnames(lhs)[(n+1):(n+k)] 

    # contRes = only continous variables (indicator variables are summed up)
    contRes <- cbind(lhs, t(out[, -1]))
    contRes <- aggregate(as.formula(paste("contRes[,-(1:(n+k))] ~", 
                                          paste(cnames, collapse = "+"))), 
                         data = contRes, 
                         na.action = na.pass,
                         FUN = sum)
    contRes <- cbind(out[, 1], t(contRes[-1, -(1:n)]))
    contNames <- c(sapply(2:(nrow(r)), function(row) 
      paste(colnames(r), r[row, 1:n], sep = "^", collapse = "*")))
    contNames <- c("time", gsub("\\*?[^\\*]+\\^0\\*?", "", contNames))
    dimnames(contRes)[[2]] <- contNames
    
    indexes <- sapply(colnames(r),
      function(s) which(stringr::str_detect(contNames, paste0("^",s, "\\^[:digit:]+$"))))

    moments <- contRes[, c(1, indexes)] %>% 
               as.data.frame() %>% 
               tidyr::gather(key = variable, value = value, -time) %>%
               dplyr::mutate(order = as.numeric(stringr::str_match(variable, "[:digit:]+$"))) %>%
               dplyr::mutate_at(.vars = "variable", 
                                .funs = stringr::str_match, pattern = "^[:alnum:]+") %>%
               tidyr::spread(key = variable, value = value) %>%
               dplyr::full_join(as.data.frame(discRes), by = "time")
    

    moments[, dnames] <- rowSums(as.matrix(moments[indicatorNames]) * 
                             sapply(states, function(i) i^moments$order))
    moments[, indicatorNames] <- NULL

    result <- structure(list(model = obj, 
                             moments = moments,
                             discRes = discRes, 
                             contRes = contRes, 
                             contInd = r,
                             maxOrder = maxOrder,
                             closure = closure,
                             centralize = centralize
                             ), class = "momApp")
    
   return(result)
})