\name{predict.bayesx}
\alias{predict.bayesx}

\title{
  Prediction from fitted BayesX MCMC model
}

\description{
  Takes a fitted \code{"bayesx"} object produced by \code{\link{bayesx}} and produces predictions
  on the scale of the \bold{linear predictor} given a new set of values for the model covariates or
  the original values used for the model fit. Note that prediction for newdata is only possible for
  \code{"MCMC"} estimated models!
}
 
\usage{
\method{predict}{bayesx}(object, newdata, model = NULL, term = NULL,
  intercept = TRUE, FUN = mean, ...)
}

\arguments{
  \item{object}{an object of class \code{"bayesx"} or \code{"bayesx.hpc"}.}
  \item{newdata}{a data frame or list containing the values of the model covariates at which
    predictions are required. Note, supplying \code{newdata} is only possible for \code{"MCMC"}
    estimated models!. If this is not provided then predictions corresponding to the original data
    are returned. If \code{newdata} is provided, then it should contain all the variables needed for 
    prediction!}
  \item{model}{for which model should predictions be calculated, either an integer or 
    a character, e.g. \code{model = "mcmc.model"}. For \code{"MCMC"} estimated models, predictions
    may be calculated by combining the parameter samples of different models. Note that the terms in
    all \code{model}s need to be exactly the same! This option is usefull when seperate models have
    been started e.g. on seperate CPUs, i.e. the user may compute predictions by combining the
    samples of all models, see also the high performance computing option (\code{hpc}) in function
    \code{\link{bayesx}}.}
  \item{term}{if not \code{NULL}, the function will search for the term the predictions should be
    calculated for, either an integer or a character, eg \code{term = "sx(x)"}. If
    \code{term = NULL} predictions for all terms will be calculated. Note that for \code{"MCMC"}
    estimated models, if more than one model \code{term} is specified, predictions will be
    calculated based on the combination of the terms specified.}
  \item{intercept}{logical. Should predictions for the intercept be included using \code{"MCMC"}
    estimated models?}
  \item{FUN}{the function that is used to calculate predictions of the samples and the data/newdata,
    e.g. could also be \code{FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975))}.}
  \item{\dots}{arguments passed to \code{\link{read.table}} if \code{newdata} is a character path
    to the data that should be used for prediction.}
}

\value{ 
  A vector or data frame with the resulting predictions. For \code{"REML"} estimated models an
  object as returned from \code{\link{fitted.bayesx}}.
}

\author{
  Nikolaus Umlauf, Thomas Kneib, Stefan Lang, Achim Zeileis.
}

\section{WARNINGS}{
  Using  geospline/geokriging model terms, predictions for MCMC models are not available at the
  moment! If geosplines are needed, it is recommended to compute the x- and y coordinates of e.g.
  the centroids of some region polygons (e.g. of a map, see function \code{\link{read.bnd}} and 
  \code{\link{plotmap}}) and use these for estimating geosplines with tensor product splines, i.e.,
  \code{bs = "te"} in function \code{\link{sx}}. In addition, out of sample prediction is not
  supported yet.
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

## get predictions for a subset of the data
nd <- dat[sample(1:n, size = 10), ]

## and only for sx(x) and w
## note that this computes the combined predictions:
## intercept + sx(x) + w
pred <- predict(b, newdata = nd, term = c("sx(x)", "w"))
pred

## predictions of multicore model
b <- bayesx(y ~ sx(x) + z + w, data = dat,
  hpc = TRUE, cores = 2)

## same as before
pred <- predict(b, newdata = nd, term = c("sx(x)", "w"))
pred
}
}

\keyword{regression}
