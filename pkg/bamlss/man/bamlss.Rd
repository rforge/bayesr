\name{bamlss}
\alias{bamlss}

\title{
  Fit Bayesian Additive Models for Location Scale and Shape (and Beyond)
}
\description{
  This is the main model fitting function of the package. Function \code{bamlss()}
  is wrapper function that parses the \code{data} and the model \code{formula}, or
  extended \code{\link{bamlss.formula}}, and the \code{\link{bamlss.family}} into a
  \code{\link{bamlss.frame}}. The \code{\link{bamlss.frame}} then holds all model
  matrices and information that is needed for setting up estimation enginges.
  The model matrices are based on \code{\link{mgcv}} infrastructures, i.e.,
  smooth terms are constructed using \code{\link{smmoth.construct}} and
  \code{\link{smoothCon}}. Therefore, all \code{\link{mgcv}} model term constructors like
  \code{\link{s}}, \code{\link{te}}, \code{\link{t2}} and \code{\link{ti}} can be used.
  Identifiablity conditions are imposed using function \code{\link{gam.side}}.

  After the \code{\link{bamlss.frame}} is set up the function applies optimizer and/or
  sampling functions to the \code{'bamlss.frame'} object. These functions can also be
  provided by the user. See the details below on how to create new engines to be used
  with function \code{bamlss()}.

  The optimum parameters and/or the samples of the parameters are then parsed to create
  summary statistics with \code{\link{summary.bamlss}}, effect plots with
  \code{\link{plot.bamlss}}, etc. See the examples.
}
\usage{
bamlss(formula, family = "gaussian", data = NULL, start = NULL,
  knots = NULL, weights = NULL, subset = NULL, offset = NULL,
  na.action = na.omit, contrasts = NULL, reference = NULL,
  transform = NULL, optimizer = NULL, sampler = NULL, samplestats = NULL,
  results = NULL, cores = NULL, sleep = NULL, combine = TRUE,
  model = TRUE, x = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{}
  \item{family}{}
  \item{data}{}
  \item{start}{}
  \item{knots}{}
  \item{weights}{}
  \item{subset}{}
  \item{offset}{}
  \item{na.action}{}
  \item{contrasts}{}
  \item{reference}{}
  \item{transform}{}
  \item{optimizer}{}
  \item{sampler}{}
  \item{samplestats}{}
  \item{results}{}
  \item{cores}{}
  \item{sleep}{}
  \item{combine}{}
  \item{model}{}
  \item{x}{}
  \item{\dots}{}
}

\details{
  Details on setting up engines \dots
}

\value{
  An object of class \code{'bamlss'}.
}

\references{
%% ~put references to the literature/web site here ~
}

\author{
  Nikolaus Umlauf, Nadja Klein, Achim Zeileis.
}

\seealso{
  \code{\link{bamlss.frame}}, \code{\link{family.bamlss}}, \code{\link{formula.bamlss}},
}

\examples{
  ## Simulated data example.
  d <- GAMart()
  b <- bamlss(num ~ s(x1) + s(x2) + s(x3), data = d)
  summary(b)
  plot(b)
}

\keyword{regression}
\keyword{MCMC}

