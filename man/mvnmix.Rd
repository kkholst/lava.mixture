\name{mvnmix}
\alias{mvnmix}
\title{Estimate mixture latent variable model}
\description{
  Estimate parameters in a mixture of latent variable models via the EM algorithm.
}
\usage{
mvnmix(data, k=2, theta, steps=500,
                 tol=1e-16, lambda=0,
                 mu=NULL,
                 silent=TRUE, extra=FALSE, ...
                 ) 
}
\arguments{
  \item{data}{\code{data.frame}}
  \item{k}{Number of mixture components}
  \item{theta}{Optional starting values}
  \item{steps}{Maximum number of iterations}
  \item{tol}{Convergence tolerance of EM algorithm}
  \item{lambda}{Added to diagonal of covariance matrix (to avoid singularities)}
  \item{mu}{Initial centres (if unspecified random centres will be
    chosen)}
  \item{silent}{Turn on/off output messages}
  \item{extra}{Extra debug information}
  \item{\dots}{Additional arguments parsed to lower-level functions}
}
\value{
  A \code{mixture} object
}
\seealso{\code{mvnmix}}
\examples{
data(faithful)
set.seed(1)
M1 <- mvnmix(faithful[,"waiting",drop=FALSE],k=2)
M2 <- mvnmix(faithful,k=2)
\donttest{
par(mfrow=c(2,1))
plot(M1,col=c("orange","blue"),ylim=c(0,0.05))
plot(M2,col=c("orange","blue"))
}
}
\author{Klaus K. Holst}
\keyword{models}
\keyword{regression}
