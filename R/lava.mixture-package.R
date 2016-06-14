

#' Estimation and simulation of mixtures of latent variable models.
#' 
#' Framwork for estimating parameters from mixtures of Latent Variable Models.
#' Plugin for the \code{lava} package.
#' 
#' 
#' @name lava.mixture-package
#' @importFrom graphics matplot plot rug par curve lines points
#' @importFrom grDevices rgb col2rgb 
#' @importFrom ucminf ucminf
#' @importFrom lava iid information score manifest parpos index<-
#' 	multigroup reindex estimate modelPar variances offdiags
#' 	Inverse Model pars covariance lvm CoefMat vars
#' @importFrom mets scoreMVN loglikMVN 
#' @importFrom stats vcov coef rmultinom runif logLik na.omit AIC
#' 	cov nlminb predict density dnorm
#' @importFrom utils tail
#' @aliases lava.mixture lava.mixture-package
#' @docType package
#' @author Klaus K. Holst Maintainer: <k.k.holst@@biostat.ku.dk>
#' @keywords package
NULL



