\name{getscript}
\alias{getscript}

\title{
  Generate an R model script  
}

\description{
  The function generates an executable \R script for obtaining summary statistics, visualization
  of model diagnostics and term effect plots of a fitted model object. 
}
 
\usage{
getscript(object, file = NULL, verbose = FALSE)
}

\arguments{
  \item{object}{an object of class \code{"bayesx"}.}
  \item{file}{optional, an output file the script is written to.}
  \item{verbose}{if set to \code{TRUE}, the script is printed to the \R console.}
}

\author{
  Thomas Kneib, Stefan Lang, Nikolaus Umlauf, Achim Zeileis.
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
b <- bayesx(y ~ s(x, bs = "ps"), data = dat)

## generate the R script
## and print it
script <- getscript(b)
script
}
}

\keyword{regression}
