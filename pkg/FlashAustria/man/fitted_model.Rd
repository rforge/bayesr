\name{flash_austria_model}
\alias{flash_austria_model}
\alias{fitted_model}

\title{A fitted BAMLSS model object}

\description{A bamlss model trained on the \code{flash_austria_train} data.
Lightning counts serve as response, the ERA5 covariates enter the model
as smooth terms using P-splines. The zero-truncated negative binomial
is employed as distribution.}

\format{
A \code{bamlss} object.
}

\references{
Umlauf, N., Klein, N., and Zeileis, A. (2018).
BAMLSS: Bayesian Additive Models for Location, Scale and Shape (and Beyond).
J. Comput. Graph. Stat., 27, 612--627.
\url{https://doi.org.10.1080/10618600.2017.1407325}
}

\examples{
data("fitted_model", package = "FlashAustria")
if (require(bamlss)) {
  plot(flash_austria_model)
}
}


\keyword{datasets}

