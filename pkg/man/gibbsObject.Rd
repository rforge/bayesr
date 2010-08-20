\name{gibbsObject}

\alias{gibbsObject}

\title{Give a detailed description of returned objects of the gibbs sampler}

\description{
 This function supplies a more detailed description of returned objects by a \code{\link{gibbs}} object.
}
 
\usage{
gibbsObject(object, ...)
}

\arguments{
  \item{object}{A returned \code{\link{gibbs}} object}
  \item{\dots}{unused}
}

\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at}}

\examples{
n<-500
x1<-runif(n,0,5)

# generate some nonlinear function
f1<-function(x){sin(x)}

# generate data 
y<-10+f1(x1)+rnorm(n,sd=1.2)

# estimate model
b<-gibbs(y~s(x1))

# show detailed description of list objects in returned gibbs object model
gibbsObject(b)

# now show names of returned gibbs object only
names(b)
}

\keyword{models}
