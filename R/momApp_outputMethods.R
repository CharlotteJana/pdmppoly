#======== todo =================================================================
#t3 in print-methode auf $moments zurückgreifen, den rest als test verwenden
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

#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_line facet_grid labs scale_color_discrete
#' @rdname momApp-methods
#' @export
plot.momApp <- function(x){
  
  data <- tidyr::gather(x$moments, key = "variable", value = "value", 
                        names(init(x$model)))
  data$order <- factor(data$order)
  
  plot <- ggplot(data = data, aes(x = time, y = value)) + 
    geom_line(aes(group = order, color = order)) +
    facet_grid(rows = ggplot2::vars(variable), scales = "free") +
    scale_color_discrete(name = "order of\nmoments") +
    labs(subtitle = paste0(descr(x$model), "\n",
                          format(x$model, short = FALSE, 
                                 slots = c("parms", "init"), collapse = "\n")),
         title = "Moment approximation",
         caption = paste("Method for moment closure:", x$closure),
         y = "")
  
  print(plot)
  invisible(plot)
  return(plot)
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