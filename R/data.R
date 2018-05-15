#======== todo =================================================================
#t1 references for every model
#t3 benaim: detailed description
#t1 PROM durch meine Promotion ersetzen
#t1 toggleSwich als vergleich mit ins example!

#========== simplePoly =====================

#' Model \code{simplePoly}
#'
#' This model is equivalent to model \code{\link{simplePdmp}}, which is given as
#' example in package \pkg{pdmpsim}. Model \code{simplePoly} is defined as a
#' \code{\link{polyPdmpModel}} object. Apart from that, it models exactly the
#' same mechanism as \code{simplePdmp}. It is included to the package for
#' demonstration purposes and is used in some unit tests and function examples.
#' @example inst/models/simplePoly.R
"simplePoly"

#========== Benaim =====================

#' Model \code{Benaim}
#' 
#' This Model was first introduced by Benaim et al in 2012 (example 1.3).
#' It is an interesting example for the importance of the switching rates.
#' @section Simulation:
#' The simulations in PROM were done with slot \code{times} set to
#' \itemize{
#' \item \code{from = 0, to = 1000, by = 0.1.}
#' }
#' The following parameter sets were simulated:
#' \itemize{
#' \item \code{b = 2, β = 1.6}
#' \item \code{b = 2, β = 1.4}
#' \item \code{b = 2, β = 1.3}
#' \item \code{b = 2, β = 0.8}
#' \item \code{b = 2, β = 0.5}
#' \item \code{b = 2, β = 0.3}
#' }
#' @example inst/models/benaim.R 
#' @format 
#'   \code{Benaim} is an object of class \code{\link[pdmpsim]{pdmpModel}}, \cr 
#'   \code{polyBenaim} is an object of class \code{\link{polyPdmpModel}}.
#' @source The model is introduced in [BenaimCo2012a] as example 1.3.
#' @name Benaim
#' @aliases benaim polyBenaim
"Benaim"

#' @rdname Benaim
"polyBenaim"

#========= genePoly 1 ===============

#' Gene regulation with constant activation
#' 
#' This PDMP models the most simple situation of gene regulation,
#' where we have one gene and a constant activation rate without
#' a further regulation mechanism. Transcription and translation
#' are considered as one step and are not modeled separately.
#' In PROM, this model is referred to as \emph{Model 1},
#' therefore it is named \code{genePdmp1} and \code{genePoly1} here.
#' @section Simulation:
#' The simulations in PROM were done with slot \code{times} set to
#' \itemize{
#' \item \code{from = 0, to = 1000, by = 0.1.}
#' }
#' The following parameter sets were simulated:
#' \itemize{
#' \item \code{κ01 = 0.01, κ10 = 0.01, α = 1, β = 0.06}
#' \item \code{κ01 = 0.01, κ10 = 0.01, α = 1, β = 0.005}
#' \item \code{κ01 = 0.01, κ10 = 0.03, α = 1, β = 0.025}
#' \item \code{κ01 = 0.03, κ10 = 0.01, α = 1, β = 0.025}
#' }
#' @example inst/models/gene1.R 
#' @format 
#'   \code{genePdmp1} is an object of class \code{\link[pdmpsim]{pdmpModel}}, \cr 
#'   \code{genePoly1} is an object of class \code{\link{polyPdmpModel}}.
#' @source The model, including most of the parameter sets, are described in
#'   [RajCo2006] and [Zeiser2009]. The parameter values do not rely on real data.
#' @name model1
#' @aliases genePoly1 gene1
"genePdmp1"

#' @rdname model1
"genePoly1"

#========= genePoly 2 ===============

#' Gene regulation with constant activation and translation
#' 
#' This PDMP models the most simple situation of gene regulation,
#' where we have one gene and a constant activation rate without
#' a further regulation mechanism. Transcription and translation
#' are modeled separately which leads to a model wit two continous
#' variables (the first (\code{ξ1}) representing the mRNA and the
#' second (\code{ξ2}) representing the protein arising from translation).
#' In PROM, this model is referred to as \emph{Model 2},
#' therefore it is named \code{genePdmp2} and \code{genePoly2} here.
#' @section Simulation:
#' The simulations in PROM were done with slot \code{times} set to
#' \itemize{
#' \item \code{from = 0, to = 1000, by = 0.1.}
#' }
#' The following parameter sets were simulated:
#' \itemize{
#' \item \code{κ01 = 0.01, κ10 = 0.01, α1 = 1, β1 = 0.06, α2 = 0.5, β2 = 0.02}
#' \item \code{κ01 = 0.01, κ10 = 0.01, α1 = 1, β1 = 0.025, α2 = 0.5, β2 = 0.02}
#' \item \code{κ01 = 0.01, κ10 = 0.03, α1 = 1, β1 = 0.025, α2 = 0.5, β2 = 0.0025}
#' }
#' @example inst/models/gene2.R 
#' @format 
#'   \code{genePdmp2} is an object of class \code{\link[pdmpsim]{pdmpModel}}, \cr 
#'   \code{genePoly2} is an object of class \code{\link{polyPdmpModel}}.
#' @source The model, including most of the parameter sets, are described in
#'   [RajCo2006] and [Zeiser2009]. The parameter values do not rely on real data.
#' @name model2
#' @aliases genePoly2 gene2
"genePdmp2"

#' @rdname model2
"genePoly2"

#========= genePoly 4 ===============

#' Gene regulation with positive feedback
#' 
#' This PDMP models the most a gene regulation mechanism where we have one gene
#' and a positive feedback loop. This means that the rate to unblock the gene
#' depends on the concentration of the gene product \code{ξ}, where a high
#' concentration leads to a higher rate and vice versa. Transcription and
#' translation are considered as one step and are not modeled separately. In
#' PROM, this model is referred to as \emph{Model 4}, therefore it is named
#' \code{genePdmp4} and \code{genePoly4} here.
#' @section Simulation:
#' The simulations in PROM were done with slot \code{times} set to
#' \itemize{
#' \item \code{from = 0, to = 1000, by = 0.1.}
#' }
#' The following parameter sets were simulated:
#' \itemize{
#' \item \code{κ01 = 0.02, κ10 = 0.02, α = 1, β = 0.2}
#' \item \code{κ01 = 0.02, κ10 = 0.02, α = 7, β = 0.2}
#' }
#' @example inst/models/gene4.R 
#' @format 
#'   \code{genePdmp4} is an object of class \code{\link[pdmpsim]{pdmpModel}}, \cr 
#'   \code{genePoly4} is an object of class \code{\link{polyPdmpModel}}.
#' @source The model, including most of the parameter sets, are described in
#'   [Zeiser2009] and [ZeiserFranzLiebscher2000]. The parameter values do not
#'   rely on real data.
#' @name model4
#' @aliases genePoly4 gene4
"genePdmp4"

#' @rdname model4
"genePoly4"

#========= genePoly 7 ===============

#' Toggle Switch with two promotors
#' 
#' This model is equivalent to model \code{\link{toggleSwitch}}, which is given
#' as example in package \pkg{pdmpsim}. Model \code{toggleSwitch} is defined as
#' a \code{\link{polyPdmpModel}} object with two discrete variables \code{dA}
#' and \code{dB}. Models \code{genePdmp7} and \code{genePoly7} model the same
#' gene regulation mechanism, describing two genes \code{A} and \code{B} that
#' mutually regulate one another. The difference is, that they are formulated
#' with only one discrete variable \code{d} that takes values 1, 2, 3, 4, where
#' \itemize{
#' \item \code{d = 1} stands for \code{dA = 0, dB = 0} (both genes are blocked),
#' \item \code{d = 2} stands for \code{dA = 1, dB = 0} (B is blocked, A is unblocked),
#' \item \code{d = 3} stands for \code{dA = 0, dB = 1} (B is unblocked, A is blocked),
#' \item \code{d = 4} stands for \code{dA = 1, dB = 1} (both genes are unblocked).
#' }
#' In PROM, the toggle switch model is referred to as \emph{Model 7}, therefore 
#' the models here are named \code{genePdmp7} and \code{genePoly7}.
#' @section Simulation:
#' The simulations in PROM were done with slot \code{times} set to
#' \itemize{
#' \item \code{from = 0, to = 1000, by = 0.1.}
#' }
#' The following parameter sets were simulated:
#' \itemize{
#' \item \code{βA = 0.02, βB = 0.02, αA = 4, αB = 4,} \cr
#' \code{κ01A = 0.05, κ10A = 0.002, κ01B = 0.05, κ10B = 0.002}
#' }
#' @example inst/models/gene7.R 
#' @format 
#'   \code{genePdmp7} is an object of class \code{\link[pdmpsim]{pdmpModel}}, \cr 
#'   \code{genePoly7} is an object of class \code{\link{polyPdmpModel}}.
#' @source The model, including most of the parameter sets, are described in
#'   [Zeiser2009]. The parameter values do not rely on real data.
#' @name model7
#' @aliases genePoly7 gene7
"genePdmp7"

#' @rdname model7
"genePoly7"

#========= genePoly 8 ===============

#' Gene regulation with positive feedback and dimerization
#' 
#' This PDMP models the most a gene regulation mechanism where we have one gene
#' and a positive feedback loop. The activator however is not the gene product
#' itself but a dimer of two molecules of the genproduct. This means that we
#' have two continous variables \code{ξ} and \code{ξd} where \code{ξ} represents
#' the gene product and \code{ξd} the concentration of the dimerized gene
#' product. Transcription and translation are considered as one step and are not
#' modeled separately. In PROM, this model is referred to as \emph{Model 8},
#' therefore it is named \code{genePdmp8} and \code{genePoly8} here.
#' @section Simulation:
#' The simulations in PROM were done with slot \code{times} set to
#' \itemize{
#' \item \code{from = 0, to = 1000, by = 0.1.}
#' }
#' The following parameter sets were simulated:
#' \itemize{
#' \item \code{κ01 = 0.02, κ10 = 0.02, α = 1, β = 0.2, γ21 = 0.1, γ12 = 0.05}
#' \item \code{κ01 = 0.02, κ10 = 0.02, α = 1, β = 0.3, γ21 = 0.1, γ12 = 0.05}
#' \item \code{κ01 = 0.02, κ10 = 0.02, α = 1, β = 0.4, γ21 = 0.1, γ12 = 0.05}
#' \item \code{κ01 = 0.02, κ10 = 0.02, α = 1, β = 0.5, γ21 = 0.1, γ12 = 0.05}
#' }
#' @example inst/models/gene8.R 
#' @format 
#'   \code{genePdmp8} is an object of class \code{\link[pdmpsim]{pdmpModel}}, \cr 
#'   \code{genePoly8} is an object of class \code{\link{polyPdmpModel}}.
#' @source The model, including most of the parameter sets, are described in
#'   [Zeiser2009]. The parameter values do not rely on real data.
#' @name model8
#' @aliases genePoly8 gene8
"genePdmp8"

#' @rdname model8
"genePoly8"
