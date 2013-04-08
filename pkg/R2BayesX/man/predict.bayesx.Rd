\name{predict.bayesx}
\alias{predict.bayesx}

\title{
  Prediction from fitted BayesX objects
}

\description{
  Takes a fitted \code{"bayesx"} object returned from \code{\link{bayesx}} and produces predictions
  by simple refitting with weights.
}
 
\usage{
\method{predict}{bayesx}(object, newdata, model = NULL, ...)
}

\arguments{
  \item{object}{an object of class \code{"bayesx"} or \code{"bayesx.hpc"}.}
  \item{newdata}{a data frame or list containing the values of the model covariates at which
    predictions are required.}
  \item{model}{for which model should predictions be calculated, either an integer or 
    a character, e.g. \code{model = "mcmc.model"}. Note that using \code{newdata} exactly one model
    must be selected within argument \code{model}!}
  \item{\dots}{arguments passed to \code{\link{fitted.bayesx}}.}
}

\value{ 
  If \code{newdata} is supplied, an object of class \code{"bayesx"} is returned, otherwise objects
  returned from \code{\link{fitted.bayesx}}. The prediction object is a simple refit of the initial
  model using the old and the \code{newdata} including \code{weights} that are set to zero for the
  \code{newdata}. Predictions can then be extracted using function \code{\link{fitted.bayesx}}.
}

\author{
  Nikolaus Umlauf, Thomas Kneib, Stefan Lang, Achim Zeileis.
}

\seealso{
  \code{\link{fitted.bayesx}}, \code{\link{bayesx}}.
}

\examples{
\dontrun{
## generate some data
set.seed(121)
n <- 500

## regressors
dat <- data.frame(x = runif(n, -3, 3), z = runif(n, 0, 1),
  w = runif(n, 0, 3))

## generate response 
dat$y <- with(dat, 1.5 + sin(x) + z -3 * w + rnorm(n, sd = 0.6))

## estimate model
b <- bayesx(y ~ sx(x) + z + w, data = dat)

## create some data for which predictions are required
nd <- data.frame(x = seq(2, 5, length = 100), z = 1, w = 0)

## prediction model from refitting with weights
p <- predict(b, newdata = nd)
plot(p)
}
}

\keyword{regression}
