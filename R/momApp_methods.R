#' Methods for objects of class \code{\link{momApp}}
#' 
#' Objects of class momApp occur as results of function \code{\link{momApp}}.
#' Their structure is described in the \code{Details} section of the
#' documentation of \code{\link{momApp}}. There are several methods to examine
#' such objects, namely \code{print}, \code{summary}, \code{plot}, \code{tail}
#' and \code{head}. Function \code{addSimulation} adds moments that come from
#' simulated data to object$moments and adds a new element \code{object$seeds}
#' containing the vector with \code{seeds} that were used for simulation. Adding
#' simulated moments makes it easy to compare the results, i.e. with
#' \code{plot}.
#' @param x object of class \code{momApp}
#' @param object object of class \code{momApp}
#' @param multSim object of class \code{\link[pdmpsim]{multSim}} that contains 
#' simulated data
#' @param ... further arguments to the default method
#' @name momApp-methods
#' @examples 
#' data(genePolyBF)
#' a <- momApp(genePolyBF, maxorder = 4, 
#'             closure = c("normal", "zero", "lognormal"), 
#'             centralize = c(TRUE, FALSE, FALSE))
#' print(a)
#' summary(a)
#' 
#' data(genePdmpBF)
#' b <- addSimulation(a, multSim(genePdmpBF, seeds = 1:30))
#' tail(b)
#' plot(b, plotorder = c(1, 4))
NULL

#' @rdname momApp-methods
#' @export
addSimulation <- function(x, multSim){
  
  if(!identical(sim(multSim$model, outSlot = FALSE, seed = 10),
                sim(x$model, outSlot = FALSE, seed = 10))){
      stop("Simulation of 'x$model' and 'multSim$model' differ, 
           they do not represent the same PDMP.")
  }
  
  msim <- NULL
  for(m in seq_len(x$maxorder)){
    msim <- dplyr::bind_rows(msim, moments(multSim, m))
  }
  msim <- cbind(method = "simulation",
                msim)
  suppressWarnings(
    x$moments <- dplyr::bind_rows(x$moments, msim)
  )
  x$moments$method <- as.factor(x$moments$method)
  x$seeds <- multSim$seeds
  return(x)
}

#' @param plotorder numerical vector specifying the orders for which
#' moments shall be plotted.
#' @param vars character vector specifying the variables
#' for which moments shall be plotted.
#' @rdname momApp-methods
#' @export
plot.momApp <- function(x, plotorder = 1, vars = names(init(x$model)), ...){
  plotdata <- reshape2::melt(x$moments, 1:3, stringsAsFactors = TRUE)
  plotdata <- subset(plotdata, order %in% plotorder)
  plotdata <- subset(plotdata, variable %in% vars)
  plotdata$variable <- as.character(plotdata$variable)

  linetypes <- c(simulation = 1,
                 `zero (raw)` = 2,
                 `zero (central)` = 3,
                 `normal (central)` = 4,
                 gamma = 5,
                 lognormal = 6,
                 `no closure` = 7)
  
  mplot <- ggplot2::ggplot(data = plotdata, ggplot2::aes(x = time, y = value)) + 
    ggplot2::geom_line(size = 1,
                       ggplot2::aes(color = method, linetype = method)) +
    ggplot2::scale_linetype_manual(values = linetypes) +
    ggplot2::scale_color_manual(values = linetypes) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                   axis.title.y = ggplot2::element_blank()) + 
    ggplot2::labs(
      title = x$model@descr,
      subtitle = format(x$model, slots = c("parms"), short = FALSE))
  
  methods <- levels(x$moments[, 1])

  if("simulation" %in% methods & length(x$seeds) > 0){
    mplot <- mplot +
      ggplot2::labs(
        caption = paste("number of simulations:", length(x$seeds)))
  }
  mplot <- mplot + 
    ggplot2::facet_wrap(variable ~ order,
              scales = "free_y", nrow = length(x$model@init),
              labeller = ggplot2::label_bquote(cols = E(.(variable)^.(order))))
  return(mplot)
}

#' @rdname momApp-methods
#' @export
print.momApp <- function(x, ...){
  cat(noquote("Model: \n"))
  cat(format(x$model, short = FALSE, collapse = "\n",
             slots = c("descr", "parms", "init")))
  
  cat("\n\nMoment approximation for moments of order > ",
      x$maxorder, " leads to \n\n", sep = "")
  
  s <- NULL
  methods <- levels(x$moments[, 1])
  
  for(m in methods){
    
    # maximal Time for which all moments are real numbers
    maxTime <-  min(vapply(seq_len(x$maxorder),
                           function(i){
                              row <- max(which(x$moments[, 1] == m &
                                               x$moments[, "order"] == i))
                              return(x$moments[row, "time"])
                           },
                    numeric(1)))
    
    s <- dplyr::bind_rows(s,
                          x$moments[which(x$moments$time == maxTime &
                                          x$moments[, 1] == m), ])
  }
  rownames(s) <- NULL
  print(s[order(s$order), ], ...)
}

#' @rdname momApp-methods
#' @export
summary.momApp <- function(object, ...){
  cat(noquote("\n$maxorder \t"), object$maxorder)
  cat(noquote("\n$closure \t"), object$closure)
  cat(noquote("\n$centralize\t"), object$centralize)
  cat(noquote("\n$model \n"))
  cat(format(object$model, short = FALSE, collapse = "\n",
             slots = c("descr", "parms", "init")))
  for(i in 1:object$maxorder){
    methods <- levels(object$moments[, 1]) 
    for(m in methods){
      cat(noquote("\n\n$moments, order = "), i, ", method = ", m, "\n")
      rows <- which(object$moments$order == i & object$moments[, 1] == m)
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
  for(closure in names(x$out)){
    cat("\n$out$", closure, "\n\n", sep = "")
    print(tail(x$out[[closure]], ...))
  }
}

#' @rdname momApp-methods
#' @importFrom utils head
#' @export
head.momApp <- function(x, ...){
  cat("$moments\n")
  print(head(x$moments, ...))
  cat("\n$out\n")
  for(closure in names(x$out)){
    cat("\n$out$", closure, "\n\n", sep = "")
    print(head(x$out[[closure]], ...))
  }
}