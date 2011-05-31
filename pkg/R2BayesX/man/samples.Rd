\name{samples}
\alias{samples}

\title{
  Extract Samples of Coefficients and Variances  
}

\description{
  Function to extract the samples generated with Markov chain Monte Carlo simulation. 
}
 
\usage{
samples(object, model = NULL, term = NULL, ...)
}

\arguments{
  \item{object}{an object of class \code{"bayesx"}.}
  \item{model}{for which model the samples should be provided, either an integer or a character, 
    e.g. \code{model = "mcmc.model"}.}
  \item{term}{if not \code{NULL}, the function will search for the term fitted values/residuals 
    specified here, either an integer or a character, eg \code{term = "s(x)"}. If set to 
    \code{term = "fixed-samples"}, the samples of the fixed effects parameters are returned}
  \item{\dots}{not used.}
}

\author{
  Thomas Kneib, Stefan Lang, Nikolaus Umlauf, Achim Zeileis.
}

\references{
Brezger, A., Lang, S. (2008): Simultaneous probability statements for Bayesian P-splines. 
  \emph{Statistical Modeling}, \bold{8}, 141--186.
}

\seealso{
  \code{\link{bayesx}}.
}

\examples{
\dontrun{
## generate some data
set.seed(111)
n <- 200

## regressor
dat <- data.frame(x = runif(n, -3, 3))

## response 
dat$y <- with(dat, 1.5 + sin(x) + rnorm(n, sd = 0.6))

## estimate model
b <- bayesx(y ~ s(x, bs = "ps"), data = dat)

## extract samples for the P-spline
sax <- samples(b, term = "s(x)")
names(sax)

## plot coefficient samples
k <- ncol(sax$coef.samples)
par(mfrow = n2mfrow(k))
for(j in 1L:k)
  plot(sax$coef.samples[,j], type = "l")

## linear effects samples
samples(b, term = "linear-samples")
}
}

\keyword{regression}
