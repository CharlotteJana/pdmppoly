#======== todo =================================================================

#' Methods for objects of class \code{\link{momApp}}
#' 
#' Objects of class momApp occur as results of function \code{\link{momApp}}.
#' Their structure is described in the \code{Details} section of the
#' documentation of \code{\link{momApp}}. There are several methods to examine
#' such objects, namely \code{print}, \code{summary}, \code{plot}, \code{tail}
#' and \code{head}. Function \code{addSimulation} allows to add moments that
#' come from simulated data to the object. This makes it easy to compare the
#' values.
#' @param x object of class \code{momApp}
#' @param object object of class \code{momApp}
#' @param ms object of class \code{\link[pdmpsim]multSim}} that contains simulated data
#' @param ... further arguments to the default method
#' @name momApp-methods
NULL

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
addSimulation <- function(x, ms){
  
}

#' @rdname momApp-methods
#' @export
print.momApp <- function(x, ...){
  cat(noquote("Model: \n"))
  cat(format(x$model, short = FALSE, collapse = "\n",
             slots = c("descr", "parms", "init")))
  
  cat("\n\nMoment approximation for moments of order > ",
      x$maxOrder, " leads to \n\n", sep = "")
  
  s <- NULL
  
  for(c in seq_along(x$closure)){
    closureName <- names(x$out[c])
    
    # maximal Time for which all moments are real numbers
    maxTime <-  min(vapply(seq_len(x$maxOrder),
                           function(i){
                              row <- max(which(x$moments[, 1] == closureName &
                                               x$moments[, "order"] == i))
                              return(x$moments[row, "time"])
                           },
                    numeric(1)))

    # cat("\n\nMoment approximation at time t = ", maxTime, " \nwith ",
    #     "moment closure method \"", x$closure[c], "\"\nfor ",
    #     ifelse(x$centralize[c], "centralized", "raw"), " moments of order > ",
    #      x$maxOrder, ": \n\n", sep = "")
    
    s <- dplyr::bind_rows(s,
                          x$moments[which(x$moments$time == maxTime &
                                          x$moments[, 1] == closureName), ])
    
    
  }
  rownames(s) <- NULL
  print(s[order(s$order), ], ...)
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
    for(c in seq_along(object$closure)){
      cat(noquote("\n\n$moments, order = "), i, ", closure = ", names(object$out)[c], "\n")
      rows <- which(object$moments$order == i & object$moments[, 1] == names(object$out)[c])
      print(summary(object$moments[rows, -(1:3)], ...))
    }
  }
}

#' @rdname momApp-methods
#' @importFrom utils tail
#' @export
tail.momApp <- function(x, ...){
  cat("$moments\n")
  print(tail(x$moments, ...))
  cat("\n$out\n")
  for(c in seq_along(x$closure)){
    cat("\n$out$", names(x$out)[c], "\n\n", sep = "")
    print(tail(x$out[[c]], ...))
  }
}

#' @rdname momApp-methods
#' @importFrom utils head
#' @export
head.momApp <- function(x, ...){
  cat("$moments\n")
  print(head(x$moments, ...))
  cat("\n$out\n")
  for(c in seq_along(x$closure)){
    cat("\n$out$", names(x$out)[c], "\n\n", sep = "")
    print(head(x$out[[c]], ...))
  }
}