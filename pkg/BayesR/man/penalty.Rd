\name{penalty}
\alias{penalty}

\title{Create a Penalty Matrix}

\description{
This function takes some design matrix and generates its related penalty matrix of an order defined by the user. It is mainly of internal
use of the functions \code{\link{ps}} and \code{\link{te}}
}

\usage{
penalty(orderpenalty, basis)
}
 
\arguments{
  \item{orderpenalty}{the order of the penalty, any order lower than the number of columns of the design matrix can be selected}
  \item{basis}{the design matrix}
}

\value{
  This function returns a penalty matrix of an order specified by the user
}

\references{ ~put references to the literature/web site here ~ }

\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}}


\examples{
# generate penalty matrix of order 2 from a p-spline design matrix
n<-500 
x<-runif(n,0,1) 

p<-ps(x)
basis<-p$basis

myPmatrix<-penalty(2,basis)
myPmatrix
}

\keyword{models} 
\keyword{smooth} 
\keyword{regression} 
