#======== todo =================================================================
#v3 contRes und discRes umbenennen
#v3 init = dirac measure?
#t1 test momApp: verschiedene closure methoden
#t1 momentClosure: documentation

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
#' @export
setGeneric("momApp",
           function(obj, maxOrder = 4, closure = "setZero")
             standardGeneric("momApp"))

#' @rdname momApp
#' @export
setMethod("momApp", signature(obj = "polyPdmpModel"), 
          function(obj, maxOrder = 4, closure = "setZero") {
  
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
    
    if(closure %in% c("normal") & !centralize){
      centralize <- TRUE
      warning("Argument 'centralize' is changed to TRUE because closure method '", closure, "'
              is only implemented for centralized moments.")
    }
    missingMoments <- momentClosure(closure, lhsMissing, lhs, n)
    
  ###### convert ode's into quoted formulas #######
    
    if(centralize){ # centralize ode's which contain missing moments
      
      # m =  highest degree that has existing ode's <- #t brauche ich das überhaupt?
      m <- max(rowSums(lhs[matchingRows[[rowsToChange]],], na.rm = TRUE)) 
      
      # schreibe die ODEs aus rowsToChange um
    }
    else{
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
    }

  #### simulate the system of odes with deSolve ######
    
    discInd <- getIndex(obj@init[dnames], states)
    state <- apply(lhs, 1, function(row) 
      if(row[n+discInd] == 1) 
        Reduce("*", obj@init[1:n]^row[1:n]) 
      else 0
    ) #initial value: dirac messure with peak in obj@init
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
                             closure = closure
                             ), class = "momApp")
    
   return(result)
})

##### help functions ####

#' A help function for \code{\link{centralMoment}}
#' 
#' This function computes the sum 
#' \deqn{\sum_{j_1 = k_1}^m ... \sum_{j_n = k_n}^m (-1)^{j-k} {p \choose j}{j \choose k}.}
#' The value of \code{centralMomentCoef} is needed in the computations in \code{\link{centralMoment}}.
#' 
#' @param p numeric vector of length n
#' @param k numeric vector of length n
#' @param m Unklar! Gilt m = p? Sind dann die Summen 0?
centralMomentCoef <- function(k, m, p){
  
  if(length(k) != length(p)) stop("length of k and p should be equal")
  if(max(k) > m) stop("m should be grater than or equal to all elements of k")
  if(min(p) < m) stop("m should be less than or equal to all elements of p")

  summationRange <- vector("list", length(k))
  for(r in seq_along(k)){
    summationRange[[r]] <- k[r]:m
  }
  summationIndexes <- do.call(expand.grid, summationRange)
  
  sum <- 0
  for(r in 1:nrow(summationIndexes)){
    index <- as.numeric(summationIndexes[r, ])
    factors <- lapply(1:length(index), function(i){
      (-1)^(index[i]-k[i])*choose(p[i], index[i])*choose(index[i], k[i])
    })
    summand <- Reduce("*", factors)
    sum <- sum + summand
  }
  return(sum)
}

#' A help function for \code{\link{momApp}}
#'
#' Let m be the maximal order for which an entry exists in lhs. 
#' Let X be a random variable and Y the corresponding centered variable, i.e. \eqn{Y = X - \mu}. 
#' This function returns the sum 
#' \deqn{E(X^p) = \sum_{k_1=0}^{p_1}...\sum_{k_n=0}^{p_n} {p \choose k} \mu^{p-k} E(Y^k),}
#' where \eqn{E(Y^k)} is replaced by an approximated value if it is not found in lhs, or
#' \deqn{E(Y^k) = \sum_{j_1=0}^{k_1}...\sum_{j_n=0}^{k_n} (-1)^{k-j} {k \choose j} \mu^{k-j} E(X^j).}
# #' \eqn{c*E(X^k)} if \eqn{m ≥ max(k)}, with
# #' \deqn{c = \sum_{j_1=k_1}^m...\sum_{j_n=k_n}^m (-1)^{j-k} {p \choose j}{j \choose k}.}
#' 
#' @inheritParams momentClosure
#' @param p numeric vector giving the order of the moment that shall be centered.
#' It should have the same length as \code{ncol(lhs)}.
#' @param k numeric vector with \code{length(k) = n}.
centralMoment <- function(lhs, n, p){
  
  m <- max(lhs[, 1:n]) # maximal order of moments for which there exist ode's
  pd <- p[(n+1):length(p)] # orders of indicator variables
  
  if(sum(p[1:n]) < m) stop("m should be less than or equal to the first n elements of p")
  # stop, falls größere Momente in p gefordert als in missingMoments enthalten sind
  # vielleicht missingMoments direkt mit momentClosure berechnen?
  
  # mu = vector with expected values
  muRowIndex <- lapply(1:n, function(i){
    muRow <- c(rep(0, n), pd)
    muRow[i] <- 1
    prodlim::row.match(muRow, lhs) #t testen, ob das geklappt hat
  })
  mu <- lapply(muRowIndex, function(i) bquote(state[.(i)])) # stimmt
  print(mu)
  
 
  k_indexes <- as.matrix(expand.grid(lapply(1:n, function(i) 0:p[i]))) # stimmt
  k_indexes_long <- cbind(k_indexes, t(replicate(nrow(k_indexes), pd)))
  names(k_indexes_long) <- names(lhs)
  s <- apply(k_indexes_long, 1, function(i) prodlim::row.match(i, lhs))
  
  lhsMissing <- k_indexes_long[is.na(s), , drop = FALSE]
  #missingMoments <- momentClosure(closure = "zero", lhs = lhs, lhsMissing = lhsMissing, n = n)
  missingMoments <- list(quote(a), quote(b), quote(c), quote(d), quote(e)) # zum testen
  
  sum <- list()
  for(k_index in 1:nrow(k_indexes)){
    k <- k_indexes[k_index, ]
    cat("k = ", k, "\n")
    if(is.na(s)[k_index]){ # if moment is missing
      index <- prodlim::row.match(c(k, pd), lhsMissing)
      pChooseK <- as.numeric(Reduce("*", lapply(1:n, function(i) choose(p[i], k[i])))) # stimmt
      muPower <- lapply(1:n, function(i) bquote(.(mu[[i]])^.(as.numeric(p[1:n]-k)[i])))
      muPower <- Reduce(function(a,b) bquote(.(a)*.(b)), muPower)
      sum[[k_index]] <- bquote(.(pChooseK)*.(muPower)*(.(missingMoments[[index]])))
    }
    else{
      inner_sum <- list()
      j_indexes <- as.matrix(expand.grid(lapply(1:n, function(i) 0:k[i])))
      
      for(j_index in 1:nrow(j_indexes)){
        j <- j_indexes[j_index, ]
        kChooseJ <- as.numeric(Reduce("*", lapply(1:n, function(i) choose(k[i], j[i]))))
        muPower <- lapply(1:n, function(i) bquote(.(mu[[i]])^.(as.numeric(k-j)[i])))
        muPower <- Reduce(function(a,b) bquote(.(a)*.(b)), muPower)
        sign <- (-1)^(sum(k)-sum(j))
        rawMomentIndex <- prodlim::row.match(c(j, pd), lhs)
        rawMoment <- bquote(state[.(rawMomentIndex)])
        inner_sum[[j_index]] <- bquote(.(sign*kChooseJ)*.(muPower)*.(rawMoment))
      }
      
      inner_sum <- Reduce(function(a,b) bquote(.(a)*.(b)), inner_sum)
      sum[[k_index]] <- bquote(.(inner_sum))
    }
  }
  sum <- Reduce(function(a,b) bquote(.(a)*.(b)), sum)
  return(sum)
  
  # if(m >= sum(k)){
  #   index <- prodlim::row.match(k, lhsMissing)
  #   pChooseK <- Reduce("*", lapply(1:n, function(i) choose(p[i], k[i])))
  #   return(bquote(.(pChooseK)*(.(missingMoments[[index]]))))
  # }
  # else{
  #   sum <- 0
  #   sumIndexes <- as.matrix(expand.grid(lapply(1:n, function(i) 0:m)))
  #   print(sumIndexes)
  #   for(r in 1:nrow(sumIndexes)){
  #     index <- sumIndexes[r, ]
  #     rawMomentIndex <- prodlim::row.match(c(index, kd), lhs)
  #     rawMoment <- bquote(state[.(rawMomentIndex)])
  #     muPower <- lapply(seq_along(mu), function(i) bquote(.(mu[[i]])^.(as.numeric(p-index)[i])))
  #     muPower <- Reduce(function(a,b) bquote(.(a)*.(b)), muPower)
  #     summand <- bquote(.(rawMoment)*.(muPower)*(.(as.numeric(closureCoefficient(index, m, p)))))
  #     sum <- bquote(.(sum) + .(summand))
  #   }
  # }
}  

#' Double factorial
#' 
#' This function returns the double factorial n!! of a natural number n.
#' It is defined as \eqn{n!! = n*(n-2)*(n-4)*...*1} if n is odd and
#' \eqn{n!! = n*(n-2)*(n-4)*...*2} if n is even.
#' 
#' @param n a non negative, natural number
#' @export
dfactorial <- function(n){
  if(n %% 1 != 0 | n < 0)
    stop("dfactorial is only defined for non negative, natural numbers")
  if(n == 0)
    return(1)
  if(n %% 2 == 1){
    odds <- seq_len(n)[as.logical(seq_len(n) %% 2)]
    return(Reduce("*", odds))
  }
  if(n %% 2 == 0){
    evens <- seq_len(n)[!as.logical(seq_len(n) %% 2)]
    return(Reduce("*", evens))
  }
}

########### moment closure ##############

#' @param n numeric: number of all continous variables
#' @param lhs  numeric matrix. Each row gives the order of a moment
#' for which an ode exists.
#' @param lhsMissing matrix. Each row gives the order of a
#' moment for which no ode exists. 
#' @param closure string specifying the closure method that is 
#' to be used. The following values are possible: \code{zero} sets
#' all moments to 0, \code{normal} sets all moments as moments of
#' a centralized normal distribution.
#' @value A list where each element is a quoted expression.
#' The i-th element of this list gives a formula for the
#' moment whose order is given in the i-th row of \code{lhsMissing}.
#' @note This function is used internally in \code{\link{momApp}}.
#' It relies on \code{lhs} and \code{lhsMissing} given in a specific
#' form: The first n columns give the moment orders of the continous
#' variables of a \code{\link{polyPdmpModel}} object, in the same order as in
#' slot \code{init}.
momentClosure <- function(closure, lhsMissing, lhs, n){
  missingMoments <- rep(list(NA), nrow(lhsMissing))
  
  # zero: all moments = 0
  if(closure == "zero")
    missingMoments <- rep(list(0), nrow(lhsMissing))
  
  # normal: moments of a centralized normal distribution
  if(closure == "normal" & n == 1){ # if we have only one continous variable
    for(i in seq_len(nrow(lhsMissing))){
      order <- as.numeric(lhsMissing[i, 1])
      sigmaRowIndex <- prodlim::row.match(c(2, lhsMissing[i, -1]), lhs) #t2 check if 2nd moment is  contained in lhs
      sigma <- bquote(sqrt(state[.(sigmaRowIndex)]))
      missingMoments[[i]] <- switch(order %% 2,
                                    bquote(.(sigma)^.(order)*.(dfactorial(order-1))),
                                    0)
    }
  }
  return(missingMoments) #t1 löschen und auf indicatorvariablen achten!!!!
  
  if(anyNA(missingMoments)){
    stop("Closure method '", closure, "' is not implemented.")
  }
  
  return(missingMoments)
}
