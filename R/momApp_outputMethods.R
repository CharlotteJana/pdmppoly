#======== todo =================================================================
#t1 Klasse momApp kommentieren

#------- output methods -----------

#' @importFrom prodlim row.match
#' @importFrom graphics par matplot legend mtext
#' @importFrom grDevices graphics.off rainbow
#' @export
plot.momApp <- function(x, ...){
  
  l <- x$degree
  k <- ncol(x$discRes) - 1
  n <- length(x$init) - 1
  
  # s = matrix with all moment combinations we are interested in
  s <- NULL
  for(i in 1:n){ 
    m <-  matrix(data = 0, nrow = l, ncol = n)
    m[,i] <- 1:l
    s <- rbind(s, m)
  }
  
  # plotRes = contRes[s]
  t <- sapply(1:nrow(s), function(i) prodlim::row.match(s[i,], x$contInd))
  plotRes <- x$contRes[,c(1, t)]
  
  # plot
  graphics.off()
  par(mar=c(6, 5, 1.5, 2), xpd=TRUE, mfrow=c(1, n), oma = c(0, 0, 2.5, 0))
  colors = rainbow(l, v = 0.85) #rainbow(nrow(s)-1, v = 0.85)
  for(i in 1:n){
    matplot(plotRes[,1], plotRes[,((i-1)*l+2):(i*l+1)], 
            lty = 1, col = colors, adj = 0, type = "l",
            xlab = "time", ylab = paste("moments of", names(x$init)[i]), ...)
    legend("bottomright", inset = c(0, -0.4), 
           fill = colors, ncol = l/2, cex = 0.75,
           legend = parse(text = colnames(plotRes[,((i-1)*l+2):(i*l+1)])))
  }
  #mtext(paste0("moment approximation for ", x$polyPdmpName), 
  #      outer = TRUE, cex = 1.2, font = 2)
  mtext(paste0("moment approximation with method \"", x$closure,"\""), 
        outer = TRUE, cex = 1.2, font = 2)
}

#' @importFrom prodlim row.match
#' @export
print.momApp <- function(x, ...){
  model <- eval(parse(text = x$polyPdmpName))
  l <- x$degree
  k <- ncol(x$discRes) - 1
  n <- length(x$init) - 1
  
  # create matrix with all moment combinations we are interested in
  s <- NULL
  for(i in 1:(n+1)){ 
    m <-  matrix(data = 0, nrow = l, ncol = n+2)
    m[,i] <- 1:l
    s <- rbind(s, m)
  }
  colnames(s) <- c(names(x$init), "moment")
  s <- as.data.frame(s)
  
  for(i in 1:(n*l)){ # fill the moments of continuous variables
    s[i, n+2] <- x$contRes[nrow(x$contRes), 
                          prodlim::row.match(as.vector(s[i,1:n]), as.matrix(x$contInd))]
  }
  for(i in 1:l){ # fill the moments of the discrete variable
    s[n*l+i, n+2] <- model@discStates[[1]]^i %*% x$discRes[nrow(x$discRes), 2:ncol(x$discRes)]
  }
  cat("Moment approximation of ", x$polyPdmpName, " up to degree ", l, " \nwith moment closure method \"", x$closure, "\" results in \n", sep = "")
  print(s[1:l,], ...)
  #write.table(format(s, justify="right"), sep = "\t", row.names=F, quote=F)
}

#' @export
summary.momApp <- function(object, ...){
  cat(noquote("$polyPdmpName \t"), object$polyPdmpName)
  cat(noquote("\n$degree \t"), object$degree)
  cat(noquote("\n$closure \t"), object$closure)
  cat(noquote("\n\n$init \n"))
  print(object$init)
  cat(noquote("\n$discRes\n"))
  print(summary(object$discRes[,-1], ...))
  cat(noquote("\n$contRes\n"))
  print(summary(object$contRes[,-1], ...))
}

#' @importFrom utils tail
#' @export
tail.momApp <- function(x, ...){
  cat(noquote("$discRes\n"))
  print(tail(x$discRes, ...))
  cat(noquote("\n$contRes\n"))
  print(tail(x$contRes, ...))
}

#' @importFrom utils head
#' @export
head.momApp <- function(x, ...){
  cat(noquote("$discRes\n"))
  print(head(x$discRes, ...))
  cat(noquote("\n$contRes\n"))
  print(head(x$contRes, ...))
}