\name{FlashAustriaModel}
\alias{FlashAustriaModel}
\alias{flash_model_ztnbinom}
\alias{flash_model_ztSICHEL}

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

\examples{\donttest{
## Visualization of fitted model
if(require("bamlss") && require("gamlss.dist") && require("gamlss.tr")) {
  data("FlashAustria", package = "FlashAustria")
  data("FlashAustriaModel", package = "FlashAustria")
  ztSICHEL <- trun(0, family = "SICHEL", local = FALSE)
  plot(flash_model_ztSICHEL)
}

## Replication code:
\dontrun{
## Required packages.
library("bamlss")
library("gamlss.dist")
library("gamlss.tr")

## Load the data.
data("FlashAustria", package = "FlashAustria")

## Generate zero truncated Sichel distribution.
ztSICHEL <- trun(0, family = "SICHEL", local = FALSE)

## Model formula, up to four parameters.
f <- list(
  counts ~ s(d2m, bs = "ps") + s(q_prof_PC1, bs = "ps") +
	s(cswc_prof_PC4, bs = "ps") + s(t_prof_PC1, bs = "ps") +
	s(v_prof_PC2, bs = "ps") + s(sqrt_cape, bs = "ps"),
         ~ s(sqrt_lsp, bs = "ps")
)

## Estimate models.
set.seed(123)
flash_model_ztnbinom <- bamlss(f,
	family = "ztnbinom", data = FlashAustriaTrain, 
	optimizer = opt_boost, cores = 3, light = TRUE, binning = TRUE, 
	maxit = 1000, thin = 3, burnin = 1000, n.iter = 2000)

set.seed(123)
flash_model_ztSICHEL <- bamlss(f,
	family = ztSICHEL, data = FlashAustriaTrain, 
	optimizer = opt_boost, cores = 3, light = TRUE, binning = TRUE, 
	maxit = 1000, thin = 3, burnin = 1000, n.iter = 2000)

}}}

\keyword{datasets}

