\name{print.gibbs}
\alias{print.gibbs}

\title{Default Printing for a Bayesian Hierarchical Generalized Structured Additive Model}

\description{
A default printing is supplied by typing only the name of the \code{\link{gibbs}} object into the command line. 
For more details see the generic function \code{\link{summary}} which will supply an elaborate summary statistic of the fitted model. 
Also type \code{names(gibbsobject)} to see all object names returned by \code{\link{gibbs}}, also see \code{\link{gibbsObject}}
}
 
\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}}

\examples{
n<-500 
x<-runif(n,-3,3)

y<-sin(x)+rnorm(n,sd=0.8)

model<-gibbs(y~ps(x))
model
}

\keyword{models} 
\keyword{smooth} 
\keyword{regression} 
