\name{FlashAustriaModel}
\alias{FlashAustriaModel}

\title{A fitted BAMLSS model object}

\description{A bamlss model trained on the \code{\link{FlashAustriaTrain}} data.
Lightning counts serve as response, the ERA5 covariates enter the model
as smooth terms using P-splines. The zero-truncated negative binomial
is employed as distribution.}

\format{
An object of class \code{"bamlss"}.
}

\references{
Umlauf, N., Klein, N., and Zeileis, A. (2018).
BAMLSS: Bayesian Additive Models for Location, Scale and Shape (and Beyond).
J. Comput. Graph. Stat., 27, 612--627.
\url{https://doi.org.10.1080/10618600.2017.1407325}
}

\examples{
data("FlashAustriaModel", package = "FlashAustria")
if(require("bamlss")) {
  plot(FlashAustriaModel)
}
}


\keyword{datasets}

