\name{FlashAustriaModel}
\alias{FlashAustriaModel}
\alias{flash_model_ztBNB}
\alias{flash_model_ztnbinom}
\alias{flash_model_ztSI}

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
  data("FlashAustriaModel", package = "FlashAustria")
  plot(flash_model_ztSI)
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
ztSI <- trun(0, family = "SI", local = FALSE)

## Generate zero truncated beta negative binomial.
ztBNB <- trun(0, family = "BNB", local = FALSE)

## Model formula, up to four parameters.
f <- list(
  counts ~ s(sqrt_cape) + s(hcc) + s(str) + s(tciw) + s(q_prof_PC1),
         ~ s(sqrt_cape) + s(hcc) + s(str) + s(tciw) + s(q_prof_PC1),
         ~ s(sqrt_cape) + s(hcc) + s(str) + s(tciw) + s(q_prof_PC1),
         ~ s(sqrt_cape) + s(hcc) + s(str) + s(tciw) + s(q_prof_PC1)
)

## Estimate models.
set.seed(123)
flash_model_ztnbinom <- bamlss(f,
  data = FlashAustriaTrain,                       
  family = "ztnbinom", binning = TRUE,            
  optimizer = opt_boost, maxit = 1000,            
  thin = 5, burnin = 1000, n.iter = 6000,         
  light = TRUE)

set.seed(123)
flash_model_ztSI <- bamlss(f,
  data = FlashAustriaTrain,
  family = ztSI, binning = TRUE,
  optimizer = opt_boost, maxit = 1000,
  thin = 5, burnin = 1000, n.iter = 6000,
  light = TRUE)

set.seed(123)
flash_model_ztBNB <- bamlss(f,
  data = FlashAustriaTrain,
  family = ztBNB, binning = TRUE,
  optimizer = opt_boost, maxit = 1000,
  thin = 5, burnin = 1000, n.iter = 6000,
  light = TRUE)
}

}}

\keyword{datasets}

