\name{n}
\alias{n}
\alias{boost.nnet.predict}

\title{Neural Networks for BAMLSS}

\description{
  This smooth constructor implements single hidden layer neural networks with sigmoidal
  activation functions.
}

\usage{
## The neural network smooth constructor.
n(..., k = 10)

## Predict function for boosted models.
boost.nnet.predict(object, ..., term = NULL, mstop = NULL)
}

\arguments{
  \item{\dots}{For function \code{n()} a formula of the type \code{~x1+x2+x3} that specifies
    the covariates that should be modeled by the neural network. For function
    \code{boost.nnet.predict()} arguments passed to function \code{\link{predict.bamlss}}.}
  \item{k}{The number of hidden nodes of the network. Note that one can set an argument
    \code{split = TRUE} to split up the neural network into, e.g., \code{nsplit = 5} parts with
    \code{k} nodes each.}
  \item{object}{A boosted \code{"bamlss"}  object.}
  \item{term}{Character, specifying the neural network model term for which predictions should be
    computed.}
  \item{mstop}{Integer, specifies the boosting iteration for which predictions should be computed.}
}

\value{
  For function \code{n()}, similar to function \code{\link[mgcv]{s}} a simple smooth
  specification object.

  For function \code{boost.nnet.predict()} predictions as returned from function
  \code{\link{predict.bamlss}}.
}


\seealso{
\code{\link{bamlss}}, \code{\link{predict.bamlss}}, \code{\link{bfit}}, \code{\link{boost}}
}

\examples{
\dontrun{## Simulate data.
set.seed(123)
d <- GAMart()

## Estimate model.
f <- num ~ n(x1) + n(x2) + n(x3) + n(~lon+lat)

b <- bamlss(f, data = d, sampler = FALSE, nu = 1)

## Plot estimated effects.
plot(b)

## Boosted version.
set.seed(111)
b <- bamlss(f, data = d, sampler = FALSE, optimizer = boost,
  nu = 0.01, stop.criterion = "BIC")

## Predict effects.
d$fx1 <- boost.nnet.predict(b,
  model = "mu", term = "n(x1)")
d$fx2 <- boost.nnet.predict(b,
  model = "mu", term = "n(x2)")
d$fx3 <- boost.nnet.predict(b,
  model = "mu", term = "n(x3)")
d$flonlat <- boost.nnet.predict(b,
  model = "mu", term = "n(~lon+lat)")

## Plot.
par(mfrow = c(2, 2))
plot2d(fx1 ~ x1, data = d)
plot2d(fx2 ~ x2, data = d)
plot2d(fx3 ~ x3, data = d)
plot3d(flonlat ~ lon + lat, data = d)
}
}

\keyword{regression}

