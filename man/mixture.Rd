\name{mixture}
\alias{mixture}
\title{Estimate mixture latent variable model}
\description{
  Estimate parameters in a mixture of latent variable models via the EM algorithm.
}
\usage{
mixture(x, data, k=length(x), control, FUN, type=c("standard","CEM","SEM"),...)
}
\arguments{
  \item{x}{List of \code{lvm} objects. If only a single \code{lvm}
    object is given, then a \code{k}-mixture of this model is fitted
    (free parameters varying between mixture components).}
  \item{k}{Number of mixture components}
  \item{control}{Optimization parameters (see details).}
  \item{FUN}{See details below.}
  \item{type}{Type of EM algorithm (standard, classification, stochastic).}
  \item{\dots}{}
}
\value{
  A \code{mixture} object
}
\author{Klaus K. Holst}
\keyword{models}
\keyword{regression}
