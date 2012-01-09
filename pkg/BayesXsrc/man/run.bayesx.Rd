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
prg <- file.path(.find.package("BayesXsrc"), "tests", "mcmc.prg")
run.bayesx(prg)
}

\keyword{regression}
