#======== todo =================================================================
#t3 in print-methode auf $moments zurückgreifen, den rest als test verwenden
#t3 in plot-methode auf $moments zurückgreifen, evtl ggplot machen
#s1 plot: Text in Legende an meine Notation anpassen

#' Methods for objects of class \code{\link{momApp}}
#' 
#' Objects of class momApp occur as results of function \code{\link{momApp}}.
#' Their structure is described in the \code{Details} section of the documentation
#' of \code{\link{momApp}}. There are several methods to examine such objects,
#' namely \code{print}, \code{summary}, \code{plot}, \code{tail} and \code{head}.
#' @param x object of class \code{momApp}
#' @param ... further arguments to the default method
#' @name momApp-methods
NULL

#------- output methods -----------

#' @importFrom prodlim row.match
#' @importFrom graphics par matplot legend mtext
#' @importFrom grDevices graphics.off rainbow
#' @rdname momApp-methods
#' @export
plot.momApp <- function(x, ...){
  
  l <- x$maxOrder
  k <- ncol(x$discRes) - 1
  n <- length(x$model@init) - length(x$model@discStates)
  
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
            xlab = "time", ylab = paste("moments of", names(x$model@init)[i]), ...)
    if(l > 1){
    legend("bottomright", inset = c(0, -0.4), 
           fill = colors, ncol = l/2, cex = 0.75,
           legend = parse(text = colnames(plotRes[,((i-1)*l+2):(i*l+1)])))
    }
  }
  mtext(paste0("moment approximation with method \"", x$closure,"\""), 
        outer = TRUE, cex = 1.2, font = 2)
}

#' @importFrom prodlim row.match
#' @rdname momApp-methods
#' @export
print.momApp <- function(x, ...){
  cat(noquote("model: \n"))
  cat(format(x$model, short = FALSE, collapse = "\n",
             slots = c("descr", "parms", "init")))
  
  l <- x$maxOrder
  k <- ncol(x$discRes) - 1
  n <- length(x$model@init) - length(x$model@discStates)
  
  # create matrix with all moment combinations we are interested in
  s <- NULL
  for(i in 1:(n+1)){ 
    m <-  matrix(data = 0, nrow = l, ncol = n+2)
    m[,i] <- 1:l
    s <- rbind(s, m)
  }
  colnames(s) <- c(names(x$model@init), "moment")
  s <- as.data.frame(s)
  
  for(i in 1:(n*l)){ # fill the moments of continuous variables
    s[i, n+2] <- x$contRes[nrow(x$contRes), 
                          prodlim::row.match(as.vector(s[i,1:n]), as.matrix(x$contInd))]
  }
  for(i in 1:l){ # fill the moments of the discrete variable
    s[n*l+i, n+2] <- x$model@discStates[[1]]^i %*% x$discRes[nrow(x$discRes), 2:ncol(x$discRes)]
  }
  cat("\n\nMoment approximation up to order ", l, " \nwith moment closure method \"", x$closure, "\" results in \n", sep = "")
  print(s, ...)
  #write.table(format(s, justify="right"), sep = "\t", row.names=F, quote=F)
}

#' @rdname momApp-methods
#' @export
summary.momApp <- function(x, ...){
  cat(noquote("\n$maxOrder \t"), x$maxOrder)
  cat(noquote("\n$closure \t"), x$closure)
  cat(noquote("\n$model \n"))
  cat(format(x$model, short = FALSE, collapse = "\n",
             slots = c("descr", "parms", "init")))
  for(i in 1:x$maxOrder){
    cat(noquote("\n\n$moments, order = "), i, "\n")
    print(summary(x$moments[which(x$moments$order == i), -(1:2)], ...))
  }
  cat(noquote("\n$discRes\n"))
  print(summary(x$discRes[,-1], ...))
  cat(noquote("\n$contRes\n"))
  print(summary(x$contRes[,-1], ...))
}

#' @rdname momApp-methods
#' @importFrom utils tail
#' @export
tail.momApp <- function(x, ...){
  cat(noquote("$moments\n"))
  print(tail(x$moments, ...))
  cat(noquote("\n$discRes\n"))
  print(tail(x$discRes, ...))
  cat(noquote("\n$contRes\n"))
  print(tail(x$contRes, ...))
}

#' @rdname momApp-methods
#' @importFrom utils head
#' @export
head.momApp <- function(x, ...){
  cat(noquote("$moments\n"))
  print(head(x$moments, ...))
  cat(noquote("$\ndiscRes\n"))
  print(head(x$discRes, ...))
  cat(noquote("\n$contRes\n"))
  print(head(x$contRes, ...))
}