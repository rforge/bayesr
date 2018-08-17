\name{boost2}
\alias{boost2}
\alias{lasso2}
\alias{bayesx2}
\alias{bboost}
\alias{predict.bboost}

\title{Some Shortcuts}

\description{
  Some simple shortcuts to model fitting engines.
}

\usage{
## BayesX.
bayesx2(...)

## Gradient boosting.
boost2(...)

## Lasso.
lasso2(...)

## Bootstrap boosting.
bboost(..., data, cores = 1, n = 2, prob = 0.6, fmstop = NULL)
\method{predict}{bboost}(object, newdata, ..., cores = 1)
}

\arguments{
  \item{data}{The data frame to be used for modeling.}
  \item{cores}{The number of cores to be used.}
  \item{n}{The number of bootstrap iterations.}
  \item{prob}{The fraction that should be used to fit the model in each bootstrap iteration.}
  \item{fmstop}{The function that should return the optimum stopping iteration. More details soon.}
  \item{object}{The \code{"bboost"} object used for prediction.}
  \item{newdata}{The data frame predictions should be made for.}
  \item{\dots}{Arguments passed to \code{\link{bamlss}} and \code{\link{predict.bamlss}}.}
}

\value{
  A \code{bamlss} object.
}


\seealso{
\code{\link{bamlss}}, \code{\link{boost}}, \code{\link{lasso}}, \code{\link{BayesX}}
}

\examples{
\dontrun{## Simulate data.
set.seed(123)
d <- GAMart()

## Estimate model.
f <- num ~ s(x1) + s(x2) + s(x3) + s(lon,lat)

## Boosted model.
b <- boost2(f, data = d)

## Plot estimated effects.
plot(b)
}
}

\keyword{regression}

