\name{cprobs}
\alias{cprobs}

\title{Calculate Contour Probabilities}

\description{
This function calculates contour probabilities of estimated functions of a \code{\link{gibbs}} object
}
 
\usage{
cprobs(x, which = 1)
}

\arguments{
  \item{x}{a \code{\link{gibbs}} object}
  \item{which}{which term should be used, e.g. the model \code{\link{formula}} looks like \code{y~x+s(z)}, if the user wants to calculate the contour probability of the 			smooth function estimated by \code{\link{s}}, then \code{which = 2} must be selected. If the model contains random effects terms modeled hierarchical, the 			\code{which} argument may be a vector of length depending on the depth of the levels specified, i.e. the model contains a third stage, then \code{which} may 			be specified as \code{which = c(2,2)}, in the case the second term of the model is a random effect and the second argument in the random effects constructor 			is modeled as another random effect level, e.g. see \code{y~-1+x+r(id1,~-1+s(z)+r(id2,~1))}. If the third stage contains an additional covariate,
		then \code{which} may be specified as \code{which = c(2,2,2)} (note that the intercept in the last stage ), 
		e.g. see the \code{\link{formula}} \code{y~-1+x+r(id1,~-1+s(z)+r(id2,~1+s(q)))}. It is also possible to calculate contour probabilities of the effects only in 			this case, simply choose \code{which = c(2,2)}, this is also valid for the second stage, \code{which = 2}
               }
}

\value{
This function returns the contour probabilty of a chosen estimated term of a \code{\link{gibbs}} object
}

\references{Brezger, A., Lang, S. (2008): Simultaneous probability statements for Bayesian P-splines. Statistical Modeling, 8, 141-186.}

\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}}

\examples{
set.seed(111)
n<-500
x1<-runif(n,0,1)
x2<-runif(n)

# generate some nonlinear function
f1<-function(x){X<-0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
                X<-X-mean(X)}

# generate data 
y<-10+f1(x1)+rnorm(n,sd=1.2)

# estimate model
b<-gibbs(y~s(x1)+s(x2))
summary(b)

# the covariate x2 should have no influence
# first check the plots
par(mfrow=c(2,1))
plot(b,which=1)
plot(b,which=2)

# now calculate contour probabilties of the terms
cprobs(b,which=1)
cprobs(b,which=2)
}

\keyword{models} 
\keyword{smooth} 
\keyword{regression} 

