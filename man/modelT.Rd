% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{modelT}
\alias{modelT}
\alias{genePdmpT}
\alias{genePolyT}
\alias{geneT}
\title{Toggle Switch with two promotors}
\format{\code{genePdmpT} is an object of class \code{\link[pdmpsim]{pdmpModel}}, \cr 
  \code{genePolyT} is an object of class \code{\link{polyPdmpModel}}.}
\source{
The model, including most of the parameter sets, are described in
  [Zeiser2009]. The parameter values do not rely on real data.
}
\usage{
genePdmpT

genePolyT
}
\description{
This model is equivalent to model \code{\link{toggleSwitch}}, which is given
as example in package \pkg{pdmpsim}. Model \code{toggleSwitch} is defined as
a \code{\link{polyPdmpModel}} object with two discrete variables \code{dA}
and \code{dB}. Models \code{genePdmpT} and \code{genePolyT} model the same
gene regulation mechanism, describing two genes \code{A} and \code{B} that
mutually regulate one another. The difference is, that they are formulated
with only one discrete variable \code{d} that takes values 1, 2, 3, 4, where
\itemize{
\item \code{d = 1} stands for \code{dA = 0, dB = 0} (both genes are blocked),
\item \code{d = 2} stands for \code{dA = 1, dB = 0} (B is blocked, A is unblocked),
\item \code{d = 3} stands for \code{dA = 0, dB = 1} (B is unblocked, A is blocked),
\item \code{d = 4} stands for \code{dA = 1, dB = 1} (both genes are unblocked).
}
In PROM, the toggle switch model is referred to as \emph{Model T}, therefore 
the models here are named \code{genePdmpT} and \code{genePolyT}.
}
\section{Simulation}{

The simulations in PROM were done with slot \code{times} set to
\itemize{
\item \code{from = 0, to = 1000, by = 0.1.}
}
The following parameter sets were simulated:
\itemize{
\item \code{bA = 0.02, bB = 0.02, aA = 4, aB = 4,} \cr
\code{k01A = 0.05, k10A = 0.002, k01B = 0.05, k10B = 0.002}
}
}

\examples{
library(spray)
#------ code to generate the pdmpModel version -----

genePdmpT <- new("pdmpModel", 
   descr = "Model T: toggleswitch with two promotors",
   parms = list(bA = 0.5, bB = 0.5, aA = 2, aB = 4, 
                 k01A = 0.5, k10A = 2, k01B = 0.3, k10B = 3),
   init = c(fA = 0.5, fB = 0.5, d = 4),
   discStates = list(d = 1:4),
   dynfunc = function(t, x, parms) {
     df <- with(as.list(c(x, parms)), 
                c(-bA * fA, -bB * fB) + switch(d, 
                                               c(0, 0), 
                                               c(aA, 0), 
                                               c(0, aB), 
                                               c(aA, aB)))
     return(c(df, 0))
   }, 
   ratefunc = function(t, x, parms) {
     return(with(as.list(c(x, parms)),
                 c(switch(d, k01B, k01B, k10B*fA, k10B*fA),
                   switch(d, k01A, k10A*fB, k01A, k10A*fB))))
   }, 
   jumpfunc = function(t, x, parms, jtype) {
     c(x[1:2], switch(jtype, 
                      switch(x[3], 3, 4, 1, 2), 
                      switch(x[3], 2, 1, 4, 3)))
   }, 
   times = c(from = 0, to = 100, by = 0.01), 
   solver = "lsodar")

#------ code to generate the polyPdmpModel version -----

library("spray")
genePolyT <- new("polyPdmpModel",
  descr = "Model T: toggleswitch with two promotors (polynomial version)",
  parms = list(bA = 0.5, bB = 0.5, aA = 2, aB = 4, 
               k01A = 0.5, k10A = 2, k01B = 0.3, k10B = 3),
  init = c(fA = 0.5, fB = 0.5, d = 4), 
  discStates = list(d = 1:4),
  dynpolys = quote(list(
    list(overall = -bA*lone(1,3), specific = list(0, aA, 0, aA)),
    list(overall = -bB*lone(2,3), specific = list(0, 0, aB, aB))
  )), 
  ratepolys = quote(list(  
    list(k01B, k01B, k10B*lone(1,3), k10B*lone(1,3)),
    list(k01A, k10A*lone(2,3), k01A, k10A*lone(2,3))
  )),
  jumpfunc = function(t, x, parms, jtype) {
    c(x[1:2], switch(jtype, 
                     switch(x[3], 3, 4, 1, 2), 
                     switch(x[3], 2, 1, 4, 3)))
  }, 
  times = c(from = 0, to = 100, by = 0.01), 
  solver = "lsodar")

#------- comparison of the models --------------

identical(sim(genePdmpT, outSlot = FALSE, seed = 10),
          sim(genePolyT, outSlot = FALSE, seed = 10))

data("toggleSwitch")
all.equal(sim(genePdmpT, outSlot = FALSE, seed = 20)[, c("fA", "fB")],
          sim(toggleSwitch, outSlot = FALSE, seed = 20)[, c("fA", "fB")],
          check.attributes = FALSE)
}
\keyword{datasets}
