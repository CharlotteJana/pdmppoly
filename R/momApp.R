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
    
  ### create all moment combinations that are needed 
  
    # r = all moment combinations that are needed
    r <- data.frame(sapply(1:n, function(i) i = 0:maxOrder))
    r <- expand.grid(r)
    if(n > 1) r <- r[which(rowSums(r) >= 0 & rowSums(r[,1:n]) <= maxOrder), ]
    dimnames(r) <- list(1:nrow(r), cnames)
    
    # create k indicator variables that indicate the state of the discrete var
    indicatorNames <- sapply(1:k, function(i) paste0(dnames, states[i]))
    indicatorMatrix <- matrix(rep(t(diag(k)), length.out = k*k*nrow(r)), 
                              ncol = k, byrow = TRUE)
   
    # t = all moment combinations times all indicator variables
    t <- as.data.frame(r[rep(seq_len(nrow(r)), each=k),])
    t <- cbind(t, indicatorMatrix)
    dimnames(t) <- list(1:nrow(t), c(cnames, indicatorNames))
    t <- as.matrix(t)
    
  ### calculate EVGenerator to every moment combination
    odeList <- lapply(1:nrow(t), function(c)
      EVGenerator(obj, m = t[c,1:n], j = states[which(t[c, (n+1):(n+k)] == 1)]))
 
  ### create system of odes  
    matchingRows <- lapply(odeList, function(ode){
      sapply(1:nrow(index(ode)), 
             function(i) prodlim::row.match(index(ode)[i,], t))
    })
    for(j in 1:length(matchingRows)){
      if(anyNA(matchingRows[[j]])){
        h <- max(rowSums(t[matchingRows[[j]],], na.rm = TRUE)) # highest degree that has existing odes
        newOde <- momentClosure(closure, odeList[[j]], h, n)
        odeList[[j]] <- newOde
        matchingRows[[j]] <- sapply(1:nrow(index(newOde)), function(i) 
          prodlim::row.match(index(newOde)[i,], t))
      }
    }
    odeSystem <- lapply(1:length(matchingRows), function(j)
        sapply(1:length(matchingRows[[j]]), function(i) {
          bquote(.(value(odeList[[j]])[i])*state[.(matchingRows[[j]][[i]])]) 
        })
    )
    odeSystem <- lapply(odeSystem, function(i) 
      Reduce(function(a,b) bquote(.(a)+.(b)), i)
    )
    
  ### simulate the system of odes with deSolve
    discInd <- getIndex(obj@init[dnames], states)
    state <- apply(t, 1, function(row) 
      if(row[n+discInd] == 1) 
        Reduce("*", obj@init[1:n]^row[1:n]) 
      else 0
    ) #initial value: dirac messure with peak in obj@init
    times <- fromtoby(obj@times)
    func <- function(t, state, parms){
      list(sapply(odeSystem, function(x) eval(x)))
    }
    out <- ode(y = state, times = times, func = func, 
               parms = obj@parms, method = obj@solver)

  ### create class 'momApp' from the result 'out' 
    # discRes = only indicator variables
    discRes <- out[, c(1:(k+1))] # Achtung: Zeitspalte kommt dazu
    colnames(discRes)[-1] <- colnames(t)[(n+1):(n+k)] 

    # contRes = only continous variables (indicator variables are summed up)
    contRes <- cbind(t, t(out[, -1]))
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

##### moment closure ####

#' @importFrom spray as.spray
#' @export
momentClosure <- function(name, ode, l, n){
  value <- value(ode)
  index <- index(ode)
  problematicRows <- which(rowSums(index) > l)
  for(i in problematicRows){
    switch(name,
           setZero = {value[i] <- 0},
           reduceDegree = {
             index[i, 1:n] <- sapply(index[i, 1:n], function(j) max(0, j-1))
             return(momentClosure(name, as.spray(index, value, addrepeats = TRUE), l, n))
             },
           stop("moment closure method \"", name, "\" is not implementet.")
    )
  }
  return(as.spray(index, value, addrepeats = TRUE))
}
