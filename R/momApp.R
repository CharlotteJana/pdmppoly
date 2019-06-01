#======== todo =================================================================
#t1 tests überarbeiten
#t1 tests so dass gamma explodiert und lognormal NaNs liefert

#' Moment approximation for polynomial PDMPs
#' 
#' This function calculates the raw moments of a polynomial PDMP. 
#' The calculation is not based on simulated data but uses the generator
#' of the PDMP to create a system of ordinary differential equations.
#' The moments can be calculated from the solution of this system of ode's. \cr \cr
#' In some cases, the system of ode's is infinetly large and can therefore
#' not be solved. A modification is necessary, meaning that all equations for 
#' moments whose order is greater than a given order will be replaced by 
#' fixed values. This process is called \emph{moment closure}. There are
#' different possibilities, which values to set. They depend on the parameter
#' \code{closure}. The following options are implemented:
#' \itemize{
#' \item \code{zero}: Set all values to zero.
#' \item \code{normal}: Set values as moments of a normal distributed variable.
#' \item \code{lognormal}: Set values as moments of a lognormal distributed variable.
#' \item \code{gamma}: Set values as moments of a gamma distributed variable.
#' }
#' Parameter \code{centralize} determines, if raw or central moments should
#' be replaced by fixed values.
#' 
#' @param obj object of class \code{\link{polyPdmpModel}}.
#' @param maxOrder integer defining the highest order of moments that are
#'   considered. Higher orders are droped and replaced by fixed values. The
#'   replacement method is specified in parameter \code{closure}.
#' @param closure character vector. Every entry defines one possibility of
#'   \emph{moment closure}, i.e. changing the system of ODEs into a closed form
#'   that is solvable. Possible values are \code{zero, normal, lognormal} or
#'   \code{gamma}.
#' @param centralize boolean vector with the same length as \code{closure}.
#'   Entry \code{centralize[i]} defines if the replacement of moments with
#'   fixed values according to \code{closure[i]} should be performed on
#'   raw or central moments.
#' @return
#' The function returns an S3-object of class \code{momApp}.
#' It contains 6 different elements:
#' \itemize{
#' \item \code{model}: the polyPdmpModel \code{obj}
#' \item \code{out}: a list. The i-th entry is the result of function
#'   \code{\link[deSolve]{ode}}, performed on the system of ode which
#'   was created with moment closure \code{closure[i]}. 
#'   Column names are added to make the result
#'   understandable.
#' \item \code{moments}: a data.frame with the resulting raw moments.
# #' of the same structure as the result of function \code{\link[pdmpsim]{moments}} in package \pkg{pdmpsim}.
#' \item \code{maxOrder}: value of parameter \code{maxOrder}
#' \item \code{closure}: value of parameter \code{closure}
#' \item \code{centralize}: value of parameter \code{centralize}
#' }
#' @examples 
#' data(genePolyBF)
#' a <- momApp(genePolyBF, maxOrder = 4)
#' plot(a)
#' print(a)
#' summary(a)
#' @name momApp
#' @aliases momentApproximation momentapproximation momapp
#' @include polypdmp_class.R polypdmp_accessors.R polypdmp_generator.R
#' @seealso \code{\link{momApp-methods}} for further analysing the result.
#' @note This method works only for PDMPs with one discrete variable.
#' @importFrom prodlim row.match
#' @importFrom simecol fromtoby
#' @importFrom deSolve ode
#' @importFrom momcalc momentList transformMoment symbolicMoments 
#' @importFrom momcalc extractCov extractMean
#' @export
setGeneric("momApp",
           function(obj, maxOrder = 4, na.rm = TRUE,
                    closure = c("zero", "zero", "normal", "lognormal", "gamma"),
                    centralize = c(TRUE, FALSE, TRUE, FALSE, FALSE))
             standardGeneric("momApp"))

#' @rdname momApp
#' @export
setMethod("momApp", signature(obj = "polyPdmpModel"), 
          function(obj, maxOrder = 4, na.rm = TRUE,
                   closure = c("zero", "zero", "normal", "lognormal", "gamma"),
                   centralize = c(TRUE, FALSE, TRUE, FALSE, FALSE)){
  
  states <- obj@discStates[[1]]
  n <- length(obj@init) - length(obj@discStates) # continuous variables
  k <- length(states)
  names <- names(obj@init)  # names of all variables
  dname <- names(obj@discStates)
  cnames <- names[!names %in% dname]
  closureName <- NULL
  
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
      closure <- "no closure"
      centralize <- NA
      closureName <- "no closure" # value of column 'closure' in the final result
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
        
        for(c in seq_along(closure)){
          if(length(missingMoments) < c)
            missingMoments[[c]] <- list()
          
          # define closureName
          if(closure[c] == "zero")
            closureName[c] <- paste0("zero (", ifelse(centralize[c], "central", "raw"), ")")
          else
            closureName[c] <- closure[c]
          if(closure[c] == "normal" & centralize[c])
            closureName[c] <- "normal (central)"
          
          if(closure[c] == "normal" & !centralize[c])
            stop("Closure method 'normal' is only implemented for centralized moments.")
          if(closure[c] %in% c("lognormal", "gamma") & centralize[c]){
            stop("Closure method '", closure[c], "' is only implemented for raw moments.")
          }
          
          # centralize ode's which contain missing moments and perform moment closure OR
          if(centralize[c]){ 
  
            mListExpanded <- momcalc::transformMoment(order = missingRow[1:n],
                                                      type = "raw",
                                                      closure = closure[c],
                                                      momentList = mList)
            
            a <- prodlim::row.match(missingRow[1:n], mListExpanded$rawMomentOrders)
            missingMoments[[c]][[i]] <- mListExpanded$rawMoments[[a]]
          }
          else{ # perform moment closure directly
            if(closure[c] == "zero"){
              cov <- NA
              mean <- NA
            }
            else{
              cov <- momcalc::extractCov(mList)
              mean <- momcalc::extractMean(mList)
            }
            missingMoments[[c]][[i]] <- momcalc::symbolicMoments(distribution = closure[c],
                                                                 missingOrders = missingRow[1:n],
                                                                 cov = cov, mean = mean)[[1]]
          }
        }
      }
    }
    
    odeSystem <- list()
    for(c in seq_along(closure)){
      odeSystem[[c]] <- rep(list(NA), nrow(lhs))
      for(j in seq_len(nrow(lhs))){
        list <- lapply(seq_along(matchingRows[[j]]), function(i){
          momentIndex <- matchingRows[[j]][[i]]
          momentIsMissing <- identical(rownames(lhsFull)[momentIndex], "missing")
          if(momentIsMissing){
            h <-  prodlim::row.match(lhsFull[momentIndex, ], lhsMissing)
            return(bquote(.(value(odeList[[j]])[i])*(.(missingMoments[[c]][[h]]))))
          }
          else{
            return(bquote(.(value(odeList[[j]])[i])*state[.(momentIndex)]))
          }
        })
        odeSystem[[c]][[j]] <- Reduce(function(a,b) bquote(.(a)+.(b)), list)
      }
    }
    
  #### simulate the system of odes with deSolve ######
    
    discInd <- getIndex(obj@init[dname], states)
    state <- apply(lhs, 1, function(row) 
        1/length(discStates(obj)[[1]])*Reduce("*", obj@init[1:n]^row[1:n]) 
    ) #initial value: dirac messure with peak in 
      # obj@init["f"] and every discrete state
    times <- fromtoby(obj@times)
    out <- list()
    for(c in seq_along(closure)){
      func <- function(lhs, state, parms){
        list(sapply(odeSystem[[c]], function(x) eval(x)))
      }
      try({
        out[[closureName[c]]] <- ode(y = state, times = times, func = func, 
                                     parms = obj@parms, method = obj@solver)
      })
    }
   
  ### create meaningful column names for out
   
   contNames <- c(sapply(2:(nrow(r)), function(row)
     paste(cnames, r[row, 1:n], sep = "^", collapse = "*")
   ))
   contNames <- gsub("\\*?[^\\*\\(]+\\^0\\*?", "", contNames) # remove ^0
   contNames <- gsub("\\^+1{1}\\*+", "\\*", contNames) # remove ^1 (part 1)
   contNames <- gsub("\\^+1{1}$", "", contNames) # remove ^1 (part 2)
   discNames <- sapply(1:k, function(i) paste0("P(", dname, "=", states[i],")"))
   outNames <- NULL
   for(i in 1:k){
     outNames <- rbind(outNames, paste0("E(", contNames, "|", 
                                       dname, "=", states[i], ")"))
   }
   outNames <- c("time", discNames, outNames)
   for(i in seq_len(length(out))){
     colnames(out[[i]]) <- outNames
   }
   
  #### create data.frame with raw moments #####
    
   moments <- data.frame()
   
   # moments of discrete variables
   colnames <- paste0("P(", dname, "=", states, ")")
   for(c in seq_along(names(out))){
     for(j in 1:maxOrder){
       values <- data.frame(method = closureName[c], order = j, time = out[[c]][, "time"])
       values[, dname] <- rowSums(out[[c]][, colnames] %*% diag(states^j))
       moments <- rbind(moments, values)
     }
   }
   
   for(name in cnames){
     moments[, name] <- NA
   }
   
   # moments of order 1
   for(i in 1:n){
     colnames <- paste0("E(", cnames[i], "|", dname, "=", states, ")")
     for(c in seq_along(names(out))){
       rows <- which(moments[, "order"] == 1 &
                     moments[, "method"] == names(out)[c] &
                     moments[, "time"] %in% out[[c]][, "time"])
       moments[rows, cnames[i]]  <- tryCatch(rowSums(out[[c]][, colnames]),
                                             error=function(cond) NA)
     }
   }
   moments$method <- factor(moments$method)
   moments$order <- factor(moments$order)
   
   # moments of order > 1
   for(i in 1:n){
     for(j in 2:maxOrder){
      colnames <- paste0("E(", cnames[i], "^", j,"|", dname, "=", states, ")")
      for(c in seq_along(names(out))){
        rows <- which(moments[, "order"] == j &
                      moments[, "method"] == names(out)[c] &
                      moments[, "time"] %in% out[[c]][, "time"])
        moments[rows, cnames[i]]  <- tryCatch(rowSums(out[[c]][, colnames]),
                                              error=function(cond) NA)
      }
     }
   }
   
   if(na.rm == TRUE)
     moments <- moments[!is.na(rowSums(moments[, -(1:3)])), ]
    
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
   