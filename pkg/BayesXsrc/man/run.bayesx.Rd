\name{run.bayesx}
\alias{run.bayesx}

\title{Run BayesX}

\description{
  Run BayesX program files from \R.
}

\usage{
run.bayesx(prg = NULL, verbose = TRUE, ...)
}

\arguments{
  \item{prg}{a file path to a \pkg{BayesX} program file. If set to \code{NULL}, \pkg{BayesX} will
    start in batch mode.}
  \item{verbose}{should output be printed to the \R console during runtime of \pkg{BayesX}.}
  \item{\dots}{further arguments to be passed to \code{\link[base]{system}}.}
}

\details{
  Function uses \code{\link[base]{system}} to run \pkg{BayesX} within an \R session.
}

\value{
  If a \code{prg} file is provided, the function returns a \code{list} containg information if
  \pkg{BayesX} was succesfully launched and how long the process was running.
}

\author{
  Daniel Adler, Thomas Kneib, Stefan Lang, Nikolaus Umlauf, Achim Zeileis.
}

\examples{
## create a temporary directory for this example
dir <- tempdir()
prg <- file.path(dir, "demo.prg")

## generate some data
set.seed(111)
n <- 200

## regressor
dat <- data.frame(x = runif(n, -3, 3))

## response
dat$y <- with(dat, 1.5 + sin(x) + rnorm(n, sd = 0.6))

## write data to dir
write.table(dat, file.path(dir, "data.raw"),
  quote = FALSE, row.names = FALSE)

## create the .prg file
writeLines("
bayesreg b
dataset d 
d.infile using data.raw
b.outfile = mcmc
b.regress y = x(psplinerw2,nrknots=20,degree=3), family=gaussian predict using d
b.getsample", prg)

## run the .prg file from R
run.bayesx(prg)
}

\keyword{regression}
