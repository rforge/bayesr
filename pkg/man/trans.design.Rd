\name{trans.design}
\alias{trans.design}

\title{Create Transformed Design Matrices and Penalties}

\description{
This function takes a design matrix \code{basis} and optionally penalty matrix \code{K} and carries out a transformation, such that the cross product of the transformed design matrices and penalties will be of diagonal form.
The transformed matrices are used in the \code{\link{gibbs}} routine to estimate a bayesian hierarchical generalized structured additive model more efficiently.
The speedup is possible since during \code{\link{gibbs}} sampling it is not necessary to draw from multivariate distributions
}

\usage{
trans.design(basis, K)
}

\arguments{
  \item{basis}{the design matrix to transform}
  \item{K}{the penalty matrix to transform}
}
 
\details{
Let \eqn{Z^{\prime}WZ = RR^{\prime}}{Z'WZ = RR'} be the Cholesky decomposition of the cross product of the design matrix. Let \eqn{QSQ^{\prime}}{QSQ'} be the singular value decomposition of 
\eqn{R^{-\prime}KR^{-1}}{R^(-1)'KR^(-1)}, then the function returns the transformed orthogonal basis \eqn{ZR^{-1}}{ZR^(-1)Q} and singular values \eqn{S}{S}. This is also known as the Demmler-Reinsch basis}

\value{
The folowing objects are returned by the function
  \item{tildeZ}{the transformed design matrix}
  \item{s}{the singular values of \eqn{R^{-\prime}KR^{-1}}{R^(-1)'KR^(-1)}}
  \item{RQ}{the matrix product \eqn{RQ}{RQ} used to retransform estimated coefficients}  
}

\references{ ~put references to the literature/web site here ~ }

\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}}

\seealso{\code{\link{chol}}, \code{\link{svd}}, \code{\link{gibbs}}}

\examples{
# generate data and some design matrix
n<-500
x<-runif(n,0,50)

p<-ps(x)

mytransform<-trans.design(p$basis,p$K)
mytransform
}

\keyword{models} 
\keyword{smooth}
\keyword{regression} 
