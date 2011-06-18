\name{plotsamples}
\alias{plotsamples}

\title{
  Plot Sampling Path(s) of Coefficient(s) and Variance(s) 
}

\description{
  This function plots the sampling paths of coefficient(s) and variance(s) stored in model term
  objects typically returned from function \code{\link{bayesx}} or \code{\link{read.bayesx.output}}.
}
 
\usage{
plotsamples(x, selected = "NA", acf = FALSE, var = FALSE, max = FALSE, ...)
}

\arguments{
  \item{x}{a vector or matrix, where each column represents a different sampling path to be 
    plotted.}
  \item{selected}{a character string containing the term name the sampling paths are plotted for.}
  \item{acf}{if set to \code{TRUE}, the autocorrelation function for each sampling path is plotted.}
  \item{var}{indicates whether coefficient or variance sampling paths are displayed and simply 
    changes the main title.}
  \item{max}{if set to \code{TRUE}, \code{plotsamples} will evaluate the maximum of each row in 
     \code{x}, i.e. the maximum of all parameters for each iteration.}
  \item{\dots}{other graphical parameters to be passed to \code{\link[graphics]{plot}} and
    \code{\link[stats]{acf}}, e.g. argument \code{lag.max} if \code{acf = TRUE}.}
}

\author{
  Thomas Kneib, Stefan Lang, Nikolaus Umlauf, Achim Zeileis.
}

\seealso{
  \code{\link{plot.bayesx}}, \code{\link{bayesx}}, \code{\link{read.bayesx.output}}.
}

\examples{
## generate some data
set.seed(111)
n <- 500

## regressors
dat <- data.frame(x = runif(n, -3, 3))

## response
dat$y <- with(dat, 1.5 + sin(x) + rnorm(n, sd = 0.6))

\dontrun{
## estimate model
b <- bayesx(y ~ s(x, bs = "ps"), data = dat)
summary(b)

## plot sampling path for
## the variance
plot(b, term = "s(x)", which = "var-samples")

## plot sampling paths for 
## coefficients
plot(b, term = "s(x)", which = "coef-samples")

## extract fitted effects
f <- fitted(b, term = "s(x)", fit.attributes = TRUE)

## now use plotsamples
cs <- attr(f, "sample")
plotsamples(cs, selected = "s(x)")

## some variations
plotsamples(cs, selected = "s(x)", acf = TRUE)
plotsamples(cs, selected = "s(x)", acf = TRUE, lag.max = 200)

## plot variance sampling path
vs <- attr(f, "variance.sample")
plotsamples(vs, selected = "s(x)", var = TRUE)
}
}

\keyword{regression}
