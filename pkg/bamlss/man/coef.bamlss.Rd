\name{coef.bamlss}
\alias{coef.bamlss}

\title{
  Extract BAMLSS Coefficients
}

\description{
  Method to extract coefficients of fitted \code{\link{bamlss}} objects, either coefficients
  returned from optimizer functions, or samples from a sampler functions.
}

\usage{
\method{coef}{bamlss}(object, model = NULL, term = NULL,
  FUN = NULL, parameters = NULL, pterms = TRUE, sterms = TRUE,
  hyper.parameters = TRUE, list = FALSE,
  full.names = TRUE, ...)
}

\arguments{
  \item{object}{An object of class \code{"bamlss"}}
  \item{model}{Character or integer. For which model should coefficients be extracted?}
  \item{term}{Character or integer. For which term should coefficients be extracted?}
  \item{FUN}{A function that is applied on the paramater samples.}
  \item{parameters}{If is set to \code{TRUE}, additionally adds estimated parameters
    returned from an optimizer function (if available).}
  \item{pterms}{Should coefficients of parameteric terms be included?}
  \item{sterms}{Should coefficients of smooths terms be included?}
  \item{hyper.parameters}{For smooth terms, should hyper parameters such as smooting variances
    be included?}
  \item{list}{Should the returned object have a list structure for each distribution parameter?}
  \item{full.names}{Should full names be assigned, indicating wether a term is paramateric "p" or
    smooth "s".}
  \item{\dots}{Arguments to be passed to \code{FUN} and function \code{\link{samples}}.}
}

\value{
  Depending on argument \code{list} and the number of distributional parameters, either a
  \code{list} or vector/matrix of model coefficients.
}

\seealso{
  \code{\link{bamlss}}.
}

\examples{
\dontrun{## Simulate data.
d <- GAMart()

## Model formula.
f <- list(
  num ~ s(x1) + s(x2) + s(x3),
  sigma ~ s(x1) + s(x2) + s(x3)
)

## Estimate model.
b <- bamlss(f, data = d)

## Extract coefficients based on MCMC samples.
coef(b)

## Now only the mean.
coef(b, FUN = mean)

## As list without the full names.
coef(b, FUN = mean, list = TRUE, full.names = FALSE)

## Coefficients only for "mu".
coef(b, model = "mu")

## And "s(x2)".
coef(b, model = "mu", term = "s(x2)")

## With optimizer parameters.
coef(b, model = "mu", term = "s(x2)", parameters = TRUE)

## Only parameteric part.
coef(b, sterms = FALSE, hyper.parameters = FALSE)

## For sigma.
coef(b, model = "sigma", sterms = FALSE,
  hyper.parameters = FALSE)
}
}

\keyword{regression}
\keyword{models}

