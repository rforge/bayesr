\name{summary.gibbs}
\alias{summary.gibbs}

\title{Summary for a Gibbs Object}

\description{
Takes a \code{\link{gibbs}} object and displays summary statistics of the sampling results
}
 
\usage{
\method{summary}{gibbs}(object, digits = 4, ...) 
}

\arguments{
  \item{object}{a fitted \code{\link{gibbs}} object}
  \item{digits}{choose the decimal places of represented numbers in the summary statistics}
  \item{\dots}{not used}
}

\details{
For any smooth functions, estimated degrees of freedom \code{edf} of the marginal functions are calculated by the trace of the hat matrix by taking either the 
mean or median (see argument \code{fitex} in the \code{\link{gibbs}} function) of the corresponding hyperparameter (smoothing parameter). Next to the \code{edf} estimated quantiles
of the hyperparameter from the MCMC sampling are displayed.
For linear effects quantiles of the estimated coefficients are provided, similar to the global variance estimation results (here \code{Sigma2}).
Additional information of the goodness of fit of the model is given by the Deviance Information Criterion \code{DIC} together with the estimated effective number of
parameters \code{pd} of the fitted model.
If the model contains random effects, the estimation results for the scale parameter are similar provided as for the global variance
}

\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}}

\seealso{\code{\link{gibbs}}}

\examples{
# generate some data
set.seed(121)
n<-500
x<-runif(n,-3,3)
z<-runif(n,0,1)
w<-runif(n,0,3)

# generate some nonlinear function
f1<-function(x){sin(x)}

# generate data 
y<-10+f1(x)+z-3*w+rnorm(n,sd=0.6)

model<-gibbs(y~ps(x)+z+w)
plot(model)

summary(model)
}

\keyword{regression}
\keyword{models}