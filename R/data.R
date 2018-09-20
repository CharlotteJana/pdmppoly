#======== todo =================================================================
#t1 references for every model
#t3 benaim: detailed description
#t1 PROM durch meine Promotion ersetzen
#t1 Raten berichtigen

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

#========= Model K ===============

#' Gene regulation with constant activation
#' 
#' This PDMP models the most simple situation of gene regulation,
#' where we have one gene and a constant activation rate without
#' a further regulation mechanism. Transcription and translation
#' are considered as one step and are not modeled separately.
#' In PROM, this model is referred to as \emph{Model K},
#' therefore it is named \code{genePdmpK} and \code{genePolyK} here.
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
#' @example inst/models/geneK.R 
#' @format 
#'   \code{genePdmpK} is an object of class \code{\link[pdmpsim]{pdmpModel}}, \cr 
#'   \code{genePolyK} is an object of class \code{\link{polyPdmpModel}}.
#' @source The model, including most of the parameter sets, are described in
#'   [RajCo2006] and [Zeiser2009]. The parameter values do not rely on real data.
#' @name modelK
#' @aliases genePolyK geneK genePdmpK
"genePdmpK"

#' @rdname modelK
"genePolyK"

#========= Model K2 ===============

#' Gene regulation with constant activation and translation
#' 
#' This PDMP models the most simple situation of gene regulation,
#' where we have one gene and a constant activation rate without
#' a further regulation mechanism. Transcription and translation
#' are modeled separately which leads to a model wit two continous
#' variables (the first (\code{ξ1}) representing the mRNA and the
#' second (\code{ξ2}) representing the protein arising from translation).
#' In PROM, this model is referred to as \emph{Model K2},
#' therefore it is named \code{genePdmpK2} and \code{genePolyK2} here.
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
#' @example inst/models/geneK2.R 
#' @format 
#'   \code{genePdmpK2} is an object of class \code{\link[pdmpsim]{pdmpModel}}, \cr 
#'   \code{genePolyK2} is an object of class \code{\link{polyPdmpModel}}.
#' @source The model, including most of the parameter sets, are described in
#'   [RajCo2006] and [Zeiser2009]. The parameter values do not rely on real data.
#' @name modelK2
#' @aliases genePolyK2 geneK2 genePdmpK2
"genePdmpK2"

#' @rdname modelK2
"genePolyK2"

#========= genePoly F ===============

#' Gene regulation with positive feedback
#' 
#' This PDMP models a gene regulation mechanism where we have one gene
#' and a positive feedback loop. This means that the rate to unblock the gene
#' depends on the concentration of the gene product \code{ξ}, where a high
#' concentration leads to a higher rate and vice versa. Transcription and
#' translation are considered as one step and are not modeled separately. In
#' PROM, this model is referred to as \emph{Model F₊}, therefore it is named
#' \code{genePdmpF} and \code{genePolyF} here.
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
#' @example inst/models/geneF.R 
#' @format 
#'   \code{genePdmpF} is an object of class \code{\link[pdmpsim]{pdmpModel}}, \cr 
#'   \code{genePolyF} is an object of class \code{\link{polyPdmpModel}}.
#' @source The model, including most of the parameter sets, are described in
#'   [Zeiser2009] and [ZeiserFranzLiebscher2000]. The parameter values do not
#'   rely on real data.
#' @name modelF
#' @aliases genePolyF genePdmpF geneF
"genePdmpF"

#' @rdname modelF
"genePolyF"

#========= genePoly KF ===============

#' Gene regulation with positive feedback and constant rates
#' 
#' This PDMP models a gene regulation mechanism similar to
#' \code{\link{genePolyF}}, where we have one gene and a positive feedback loop.
#' The rate to unblock the gene depends on the concentration of the gene product
#' \code{ξ}, but it is never zero because there is an additional rate that is
#' independet of \code{ξ}. Transcription and translation are considered as one
#' step and are not modeled separately. In PROM, this model is referred to as
#' \emph{Model F₊}, therefore it is named \code{genePdmpKF} and
#' \code{genePolyKF} here.
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
#' @example inst/models/geneKF.R 
#' @format 
#'   \code{genePdmpKF} is an object of class \code{\link[pdmpsim]{pdmpModel}},\cr 
#'   \code{genePolyKF} is an object of class \code{\link{polyPdmpModel}}.
#' @source The parameter values do not rely on real data.
#' @name modelKF
#' @aliases genePolyKF genePdmpKF geneKF
"genePdmpKF"

#' @rdname modelKF
"genePolyKF"

#========= genePoly BF ===============

#' Gene regulation with positive feedback with basal transcription
#' 
#' This PDMP models a gene regulation mechanism similar to
#' \code{\link{genePolyF}}, where we have one gene and a positive feedback loop.
#' The difference is that in both discrete states transcription takes place, but
#' with different rates \eqn{α_0, α_1}, where \eqn{α_0 < α_1}. Transcription and
#' translation are considered as one step and are not modeled separately. In
#' PROM, this model is referred to as \emph{Model BF₊}, therefore it is named
#' \code{genePdmpBF} and \code{genePolyBF} here.
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
#' @example inst/models/geneBF.R 
#' @format 
#'   \code{genePdmpBF} is an object of class \code{\link[pdmpsim]{pdmpModel}},\cr 
#'   \code{genePolyBF} is an object of class \code{\link{polyPdmpModel}}.
#' @source The parameter values do not rely on real data.
#' @name modelBF
#' @aliases genePolyBF genePdmpBF geneBF
"genePdmpBF"

#' @rdname modelBF
"genePolyBF"

#========= Model T ===============

#' Toggle Switch with two promotors
#' 
#' This model is equivalent to model \code{\link{toggleSwitch}}, which is given
#' as example in package \pkg{pdmpsim}. Model \code{toggleSwitch} is defined as
#' a \code{\link{polyPdmpModel}} object with two discrete variables \code{dA}
#' and \code{dB}. Models \code{genePdmpT} and \code{genePolyT} model the same
#' gene regulation mechanism, describing two genes \code{A} and \code{B} that
#' mutually regulate one another. The difference is, that they are formulated
#' with only one discrete variable \code{d} that takes values 1, 2, 3, 4, where
#' \itemize{
#' \item \code{d = 1} stands for \code{dA = 0, dB = 0} (both genes are blocked),
#' \item \code{d = 2} stands for \code{dA = 1, dB = 0} (B is blocked, A is unblocked),
#' \item \code{d = 3} stands for \code{dA = 0, dB = 1} (B is unblocked, A is blocked),
#' \item \code{d = 4} stands for \code{dA = 1, dB = 1} (both genes are unblocked).
#' }
#' In PROM, the toggle switch model is referred to as \emph{Model T}, therefore 
#' the models here are named \code{genePdmpT} and \code{genePolyT}.
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
#' @example inst/models/geneT.R 
#' @format 
#'   \code{genePdmpT} is an object of class \code{\link[pdmpsim]{pdmpModel}}, \cr 
#'   \code{genePolyT} is an object of class \code{\link{polyPdmpModel}}.
#' @source The model, including most of the parameter sets, are described in
#'   [Zeiser2009]. The parameter values do not rely on real data.
#' @name modelT
#' @aliases genePolyT geneT genePdmpT
"genePdmpT"

#' @rdname modelT
"genePolyT"

#========= genePoly DF ===============

#' Gene regulation with positive feedback and dimerization
#' 
#' This PDMP models the most a gene regulation mechanism where we have one gene
#' and a positive feedback loop. The activator however is not the gene product
#' itself but a dimer of two molecules of the genproduct. This means that we
#' have two continous variables \code{ξ} and \code{ξd} where \code{ξ} represents
#' the gene product and \code{ξd} the concentration of the dimerized gene
#' product. Transcription and translation are considered as one step and are not
#' modeled separately. In PROM, this model is referred to as \emph{Model DF},
#' therefore it is named \code{genePdmpDF} and \code{genePolyDF} here.
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
#' @example inst/models/geneDF.R 
#' @format 
#'   \code{genePdmpDF} is an object of class \code{\link[pdmpsim]{pdmpModel}}, \cr 
#'   \code{genePolyDF} is an object of class \code{\link{polyPdmpModel}}.
#' @source The model, including most of the parameter sets, are described in
#'   [Zeiser2009]. The parameter values do not rely on real data.
#' @name modelDF
#' @aliases genePolyDF genePdmpDF geneDF
"genePdmpDF"

#' @rdname modelDF
"genePolyDF"
