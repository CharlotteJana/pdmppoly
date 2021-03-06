% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{modelK}
\alias{modelK}
\alias{genePdmpK}
\alias{genePolyK}
\alias{geneK}
\title{Gene regulation with constant activation}
\format{\code{genePdmpK} is an object of class \code{\link[pdmpsim]{pdmpModel}}, \cr 
  \code{genePolyK} is an object of class \code{\link{polyPdmpModel}}.}
\source{
The model, including most of the parameter sets, are described in
  [RajCo2006] and [Zeiser2009]. The parameter values do not rely on real data.
}
\usage{
genePdmpK

genePolyK
}
\description{
This PDMP models the most simple situation of gene regulation,
where we have one gene and a constant activation rate without
a further regulation mechanism. Transcription and translation
are considered as one step and are not modeled separately.
In PROM, this model is referred to as \emph{Model K},
therefore it is named \code{genePdmpK} and \code{genePolyK} here.
}
\section{Simulation}{

The simulations in PROM were done with slot \code{times} set to
\itemize{
\item \code{from = 0, to = 1000, by = 0.1.}
}
The following parameter sets were simulated:
\itemize{
\item \code{k01 = 0.01, k10 = 0.01, a = 1, b = 0.06}
\item \code{k01 = 0.01, k10 = 0.01, a = 1, b = 0.005}
\item \code{k01 = 0.01, k10 = 0.03, a = 1, b = 0.025}
\item \code{k01 = 0.03, k10 = 0.01, a = 1, b = 0.025}
}
}

\examples{
library(spray)
#------ code to generate the pdmpModel version -----

genePdmpK <- new("pdmpModel",
   descr = "Model K: constant activation",
   parms = list(b = 0.005, a = 1, k10 = 0.01, k01 = 0.01),
   init = c(f = 0, d = 1),
   discStates = list(d = 0:1),
   dynfunc = function(t, x, parms) {
     df <- with(as.list(c(x, parms)), a*d - b*f)
     return(c(df, 0))
   },
   ratefunc = function(t, x, parms) {
     return(with(as.list(c(x, parms)), switch(d + 1, k01, k10)))
   },
   jumpfunc = function(t, x, parms, jtype) {
     c(x[1], 1 - x[2])
   },
   times = c(from = 0, to = 100, by = 0.1),
   solver = "lsodar")

#------ code to generate the polyPdmpModel version -----

genePolyK <- new("polyPdmpModel",
     descr = "Model K: constant activation (polynomial version)",
     parms = list(b = 0.005, a = 1, k10 = 0.01, k01 = 0.01),
     init = c(f = 0, d = 1),
     discStates = list(d = 0:1),
     dynpolys = quote(list(
       list(overall = linear(c(-b, a)))
     )),
     ratepolys = quote(list(
       list(k01, k10)
     )),
     jumpfunc = function(t, x, parms, jtype) {
       c(x[1], 1 - x[2])
     },
     times = c(from = 0, to = 100, by = 0.1),
     solver = "lsodar")

#------- comparison of the models --------------

identical(sim(genePdmpK, outSlot = FALSE, seed = 20),
          sim(genePolyK, outSlot = FALSE, seed = 20))
}
\keyword{datasets}
