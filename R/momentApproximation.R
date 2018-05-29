#======== todo =================================================================
# Klasse momApp kommentieren
# warum werden bei str(...) bei contInd so komische Sachen angezeigt?
# contRes und discRes umbenennen
# init rausnehmen?
# use messages (auch bei "all.equal was used")
# Aufruf mit l=1 gibt Fehler

#' Moment approximation for polynomial PDMPs
#' 
#' @param obj object of class \code{\link{polyPdmpModel}}.
#' @param l integer defining the highest degree of moments that are considered.
#'   Higher degrees are droped and replaced by other values. The replacement
#'   method is specified in parameter \code{closure}.
#' @param closure string defining the method that does the moment closure, i. e.
#'   that changes the system of ODEs into a closed form that is solvable.
#'   Possible values are \code{setZero} (the default) and reduceDegree.
#' @param plot boolean value indicating if the results shall be plotted.
#'   Defaults to TRUE.
#' @name momentApprox
#' @aliases momentApproximation momentapproximation momentapprox
#' @include polypdmp_class.R polypdmp_accessors.R polypdmp_generator.R
#' @seealso \code{link{momentClosure}} for the internal method that performs 
#'   the moment closure.
#' @return an object of class \code{\link{momApp}}.
#' @export
setGeneric("momentApprox",
           function(obj, l = 4, closure = "setZero", plot = TRUE, ...)
             standardGeneric("momentApprox"))

#' @rdname momentApprox
#' @export
setMethod("momentApprox", signature(obj = "polyPdmpModel"), 
          function(obj, l = 4, closure = "setZero", plot = TRUE) {
  
  n <- length(obj@init) - length(obj@discStates) # continuous variables
  k <- length(obj@discStates[[1]])
  names <- names(obj@init)  # names of all variables
    
  ### create all moment combinations that are needed 
  
    # r = all moment combinations that are needed
    r <- data.frame(sapply(1:n, function(i) i = 0:l))
    r <- expand.grid(r)
    if(n > 1) r <- r[which(rowSums(r) >= 0 & rowSums(r[,1:n]) <= l), ]
    dimnames(r) <- list(1:nrow(r), names[-(n+1)])
    
    # create k indicator variables that indicate the state of the discrete variable 
    indicatorNames <- sapply(1:k, function(i) paste(names[n+1], obj@discStates[[1]][i], sep = ""))
    indicatorMatrix <- matrix(rep(t(diag(k)), length.out = 2*k*nrow(r)), ncol = k, byrow = TRUE)
   
    # t = all moment combinations times all indicator variables
    t <- as.data.frame(r[rep(seq_len(nrow(r)), each=k),])
    t <- cbind(t, indicatorMatrix)
    dimnames(t) <- list(1:nrow(t), c(names[1:n], indicatorNames))
    t <- as.matrix(t)
    
  ### calculate EVGenerator to every moment combination
    odeList <- lapply(1:nrow(t), function(c)
      EVGenerator(obj, m = t[c, 1:n], i = which(t[c, (n+1):(n+k)] == 1)))
 
  ### create system of odes  
    matchingRows <- lapply(odeList, function(ode){
      sapply(1:nrow(index(ode)), function(i) row.match(index(ode)[i,], t))
    })
    for(j in 1:length(matchingRows)){
      if(anyNA(matchingRows[[j]])){
        h <- max(rowSums(t[matchingRows[[j]],], na.rm = TRUE)) # highest degree that has existing odes
        newOde <- momentClosure(closure, odeList[[j]], h, n)
        odeList[[j]] <- newOde
        matchingRows[[j]] <- sapply(1:nrow(index(newOde)), function(i) row.match(index(newOde)[i,], t))
      }
    }
    odeSystem <- lapply(1:length(matchingRows), function(j)
        sapply(1:length(matchingRows[[j]]), function(i) {
          bquote(.(value(odeList[[j]])[i])*state[.(matchingRows[[j]][[i]])]) 
        })
    )
    odeSystem <- lapply(odeSystem, function(i) Reduce(function(a,b) bquote(.(a)+.(b)), i))
    
  ### simulate the system of odes with deSolve
    discInd <- getIndex(obj@init[length(obj@init)], obj@discStates[[1]])
    state <- apply(t, 1, function(row) if(row[n+discInd] == 1) Reduce("*", obj@init[1:n]^row[1:n]) else 0) #initial value: dirac messure with peak in obj@init
    times <- fromtoby(obj@times)
    func <- function(t, state, parms){
      list(sapply(odeSystem, function(x) eval(x)))
    }
    out <- ode(y = state, times = times, func = func, parms = obj@parms, method = obj@solver)
     
  ### create class 'momApp' from the result 'out' 
  
    # discRes = only indicator variables
    discRes <- out[, c(1:(k+1))] # Achtung: Zeitspalte kommt dazu
    colnames(discRes)[-1] <- colnames(t)[(n+1):(n+k)] 
    
    # contRes = only continous variables (indicator variables are summed up)
    contRes <- cbind(t, t(out[, -1]))
    contRes <- aggregate(as.formula(paste("contRes[,-(1:(n+k))] ~", paste(colnames(contRes)[1:n], collapse = "+"))), data=contRes, FUN=sum)
    #contInd <- as.matrix(contRes[, 1:n], dimnames = list(NULL, names[1:n]))
    contRes <- cbind(out[, 1], t(contRes[-1, -(1:n)]))
    contNames <- c(sapply(2:(nrow(r)), function(row) paste(colnames(r), r[row, 1:n], sep = "^", collapse = "*")))
    contNames <- c("time", gsub("\\*?[^\\*]+\\^0\\*?", "", contNames))
    dimnames(contRes)[[2]] <- contNames
  
    result <- structure(list(polyPdmpName = deparse(eval(substitute(substitute(obj)), parent.frame())), 
                             discRes = discRes, 
                             contRes = contRes, 
                             contInd = r,
                             init = obj@init,
                             degree = l,
                             closure = closure
                             ), class = "momApp")
 
  ### return 
   if(plot) plot(result)
   return(result)
})

##### moment closure ####

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
           {stop("moment closure method \"", name, "\" is not implementet.")}
    )
  }
  return(as.spray(index, value, addrepeats = TRUE))
}
#### output methods ####

plot.momApp <- function(object){
  
  l <- object$degree
  k <- ncol(object$discRes) - 1
  n <- length(object$init) - 1
  
  # s = matrix with all moment combinations we are interested in
  s <- NULL
  for(i in 1:n){ 
    m <-  matrix(data = 0, nrow = l, ncol = n)
    m[,i] <- 1:l
    s <- rbind(s, m)
  }
  
  # plotRes = contRes[s]
  t <- sapply(1:nrow(s), function(i) row.match(s[i,], object$contInd))
  plotRes <- object$contRes[,c(1, t)]
  
  # plot
  graphics.off()
  par(mar=c(6, 5, 1.5, 2), xpd=TRUE, mfrow=c(1, n), oma = c(0, 0, 2.5, 0))
  colors = rainbow(l, v = 0.85) #rainbow(nrow(s)-1, v = 0.85)
  for(i in 1:n){
    matplot(plotRes[,1], plotRes[,((i-1)*l+2):(i*l+1)], lty = 1, col = colors, adj = 0,
            xlab = "time", ylab = paste("moments of", names(object$init)[i]), type = "l")
    legend("bottomright", inset = c(0, -0.4), fill = colors, ncol = l/2, cex = 0.75,
           legend = parse(text = colnames(plotRes[,((i-1)*l+2):(i*l+1)])))
  }
  #mtext(paste("moment approximation for ", object$polyPdmpName, sep = ""), outer = TRUE, cex = 1.2, font = 2)
  mtext(paste("moment approximation with method \"", object$closure,"\"", sep = ""), outer = TRUE, cex = 1.2, font = 2)
  
}
print.momApp <- function(object){
  model <- eval(parse(text = object$polyPdmpName))
  l <- object$degree
  k <- ncol(object$discRes) - 1
  n <- length(object$init) - 1
  
  # create matrix with all moment combinations we are interested in
  s <- NULL
  for(i in 1:(n+1)){ 
    m <-  matrix(data = 0, nrow = l, ncol = n+2)
    m[,i] <- 1:l
    s <- rbind(s, m)
  }
  colnames(s) <- c(names(object$init), "moment")
  s <- as.data.frame(s)
  
  for(i in 1:(n*l)){ # fill the moments of continuous variables
    s[i, n+2] <- object$contRes[nrow(object$contRes), row.match(as.vector(s[i,1:n]), as.matrix(object$contInd))]
  }
  for(i in 1:l){ # fill the moments of the discrete variable
    s[n*l+i, n+2] <- model@discDomain^i %*% object$discRes[nrow(object$discRes), 2:ncol(object$discRes)]
  }
  cat("Moment approximation of ", object$polyPdmpName, " up to degree ", l, " \nwith moment closure method \"", object$closure, "\" results in \n", sep = "")
  print(s[1:l,])
  #write.table(format(s, justify="right"), sep = "\t", row.names=F, quote=F)
}
summary.momApp <- function(object){
  cat(noquote("$polyPdmpName \t"), object$polyPdmpName)
  cat(noquote("\n$degree \t"), object$degree)
  cat(noquote("\n$closure \t"), object$closure)
  cat(noquote("\n\n$init \n"))
  print(object$init)
  cat(noquote("\n$discRes\n"))
  print(summary(object$discRes[,-1]))
  cat(noquote("\n$contRes\n"))
  print(summary(object$contRes[,-1]))
}
tail.momApp <- function(object){
  cat(noquote("$discRes\n"))
  print(tail(object$discRes))
  cat(noquote("\n$contRes\n"))
  print(tail(object$contRes))
}
head.momApp <- function(object){
  cat(noquote("$discRes\n"))
  print(head(object$discRes))
  cat(noquote("\n$contRes\n"))
  print(head(object$contRes))
}