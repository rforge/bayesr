\name{FlashAustriaModel}
\alias{FlashAustriaModel}

\title{Fitted Distributional Regression Model Object (bamlss)}

\description{A zero-truncated negative binomial \code{\link[bamlss]{bamlss}}
model trained on the \code{\link{FlashAustriaTrain}} data using
gradient boosting with subsequent MCMC sampling. The response variable
are the lightning counts and the regression terms are P-splines
based on ERA5 covariates.}

\format{
An object of class \code{\link[bamlss]{bamlss}}.
}

\references{
Umlauf, N., Klein, N., and Zeileis, A. (2018).
BAMLSS: Bayesian Additive Models for Location, Scale and Shape (and Beyond).
Journal of Computational and Graphical Statistics, 27(3), 612--627.
\doi{10.1080/10618600.2017.1407325}
}

\examples{
## Visualization of fitted model
if(require("bamlss")) {
data("FlashAustriaModel", package = "FlashAustria")
plot(FlashAustriaModel)
}

## Replication code:
\dontrun{
## Learning data
data("FlashAustria", package = "FlashAustria")

## Model formula
f <- list(
  counts ~ s(d2m, bs = "ps") + s(q_prof_PC1, bs = "ps") +
    s(cswc_prof_PC4, bs = "ps") + s(t_prof_PC1, bs = "ps") +
    s(v_prof_PC2, bs = "ps") + s(sqrt_cape, bs = "ps"), 
  theta ~ s(sqrt_lsp, bs = "ps")
)

## Model fitting using boosting with subsequent MCMC sampling
set.seed(111)
FlashAustriaModel <- bamlss(f,
  family = "ztnbinom", data = FlashAustriaTrain,
  optimizer = boost, maxit = 1000,
  thin = 5, burnin = 1000, n.iter = 6000)
}
}

\keyword{datasets}

