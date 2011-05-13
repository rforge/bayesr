\name{cprob}
\alias{cprob}

\title{
  Extract Contour Probabilities  
}

\description{
  Function to extract estimated contour probabilities of a particular effect estimated with 
  P-splines using Markov chain Monte Carlo estimation techniques. 
}
 
\usage{
cprob(object, model = NULL, term = NULL, ...)
}

\arguments{
  \item{object}{an object of class \code{"bayesx"}.}
  \item{model}{for which model the fitted values/residuals should be provided, either an integer or 
    a character, e.g. \code{model = "mcmc.model"}.}
  \item{term}{if not \code{NULL}, the function will search for the term fitted values/residuals 
    specified here, either an integer or a character, eg \code{term = "s(x)"}.}
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
n <- 500

## regressor
dat <- data.frame(x = runif(n, -3, 3))

## response 
dat$y <- with(dat, 1.5 + sin(x) + rnorm(n, sd = 0.6))

## estimate model
b <- bayesx(y ~ s(x, bs = "ps", xt = list(contourprob = 4)), data = dat)

## extract contour probabilities
cprob(b, term = "s(x)")
}
}

\keyword{regression}
