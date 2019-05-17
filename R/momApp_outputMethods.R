#======== todo =================================================================

#' Methods for objects of class \code{\link{momApp}}
#' 
#' Objects of class momApp occur as results of function \code{\link{momApp}}.
#' Their structure is described in the \code{Details} section of the
#' documentation of \code{\link{momApp}}. There are several methods to examine
#' such objects, namely \code{print}, \code{summary}, \code{plot}, \code{tail}
#' and \code{head}.
#' @param x object of class \code{momApp}
#' @param object object of class \code{momApp}
#' @param ... further arguments to the default method
#' @name momApp-methods
NULL

#------- output methods -----------

#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_line facet_grid labs scale_color_discrete
#' @rdname momApp-methods
#' @method plot momApp
#' @export
plot.momApp <- function(x, ...){
  
  # to avoid the R CMD Check NOTE 'no visible binding for global variable ...'
  time <- variable <- init <- descr <- NULL
  
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

#' @rdname momApp-methods
#' @export
print.momApp <- function(x, ...){
  cat(noquote("model: \n"))
  cat(format(x$model, short = FALSE, collapse = "\n",
             slots = c("descr", "parms", "init")))
  
  # maximal Time for which all moments are real numbers
  maxTime <- min(vapply(seq_len(maxOrder),
                    function(i) {
                      h <- x$moments[which(x$moments[, "order"] == i), ]
                      pos <- Position(is.na, rowSums(h))
                      ifelse(is.na(pos),
                             times(x$model)["to"],
                             h[, "time"][pos - 1])
                    },
                    numeric(1)))
  
  cat("\n\nMoment approximation at time t = ", maxTime, " \nwith ",
      "moment closure method \"", x$closure, "\"\nfor ",
      ifelse(x$centralize, "centralized", "raw"), " moments of order > ",
      x$maxOrder, ": \n\n", sep = "")
  
  s <- x$moments[which(x$moments$time == maxTime), ]
  rownames(s) <- NULL
  s$time <- NULL
  print(s, ...)
}

#' @rdname momApp-methods
#' @export
summary.momApp <- function(object, ...){
  cat(noquote("\n$maxOrder \t"), object$maxOrder)
  cat(noquote("\n$closure \t"), object$closure)
  cat(noquote("\n$centralize\t"), object$centralize)
  cat(noquote("\n$model \n"))
  cat(format(object$model, short = FALSE, collapse = "\n",
             slots = c("descr", "parms", "init")))
  for(i in 1:object$maxOrder){
    cat(noquote("\n\n$moments, order = "), i, "\n")
    print(summary(object$moments[which(object$moments$order == i), -(1:2)], ...))
  }
  cat(noquote("\n$discRes\n"))
  print(summary(object$discRes[,-1], ...))
  cat(noquote("\n$contRes\n"))
  print(summary(object$contRes[,-1], ...))
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