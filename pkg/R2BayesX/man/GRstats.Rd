\name{GRstats}
\alias{GRstats}

\title{
  Compute Gelman and Rubin's convergence diagnostics from multicore BayesX models.
}

\description{
  This function takes a fitted \code{\link{bayesx}} object estimated using the high performance
  option and computes the Gelman and Rubin's convergence diagnostic of the model parameters using
  function \code{\link[coda]{gelman.diag}} provided in package \pkg{coda}. Note that
  high performance computing is not available on Windows systems, see function
  \code{\link[parallel]{mclapply}}.
}
 
\usage{
GRstats(object, term = NULL, combine = TRUE, ...)
}

\arguments{
  \item{object}{an object of class \code{"bayesx.hpc"}, returned from the model fitting function
    \code{\link{bayesx}} using the high performance computing option (\code{hpc}).}
  \item{term}{character or integer. The term for which the diagnostics should be computed, 
    see also function \code{\link{samples}}.}
  \item{combine}{logical. Should the diagnostics be computed for every \code{term} separately?}
  \item{\dots}{arguments passed to function \code{\link[coda]{gelman.diag}}.}
}

\author{
  Nikolaus Umlauf, Thomas Kneib, Stefan Lang, Achim Zeileis.
}

\seealso{
  \code{\link{bayesx}}, \code{\link[coda]{gelman.diag}}, \code{\link{samples}}.
}

\examples{
\dontrun{
## generate some data
set.seed(111)
n <- 500

## regressors
dat <- data.frame(x = runif(n, -3, 3), z = runif(n, -3, 3),
   w = runif(n, 0, 6), fac = factor(rep(1:10, n/10)))

## response
dat$y <- with(dat, 1.5 + sin(x) + cos(z) * sin(w) +
   c(2.67, 5, 6, 3, 4, 2, 6, 7, 9, 7.5)[fac] + rnorm(n, sd = 0.6))

## estimate model
b <- bayesx(y ~ sx(x) + sx(z, w, bs = "te") + fac,
   data = dat, method = "MCMC", hpc = TRUE, cores = 2)

## obtain Gelman and Rubin's convergence diagnostics
GRstats(b, term = c("sx(x)", "sx(z,w)"))
GRstats(b, term = c("linear-samples", "var-samples"))

## of all parameters
GRstats(b, term = c("sx(x)", "sx(z,w)",
  "linear-samples", "var-samples"))
}
}

\keyword{regression}
