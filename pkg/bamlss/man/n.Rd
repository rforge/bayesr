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
b <- bamlss(f, data = d, sampler = FALSE, optimizer = boost, nu = 0.01, maxit = 1000)

## Predict.
nd <- unique(d[, c("lon", "lat")])
for(i in c(10, 20, 100, 1000)) {
  nd[[paste("p", i, sep = "")]] <- boost.nnet.predict(b, newdata = nd,
    model = "mu", term = "n(~lon+lat)", mstop = i)
}

## Plot.
par(mfrow = c(2, 2))
plot3d(p10 ~ lon + lat, data = nd)
plot3d(p20 ~ lon + lat, data = nd)
plot3d(p100 ~ lon + lat, data = nd)
plot3d(p1000 ~ lon + lat, data = nd)
}
}

\keyword{regression}

