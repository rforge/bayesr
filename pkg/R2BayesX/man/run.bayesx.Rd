\name{run.bayesx}
\alias{run.bayesx}

\title{Run BayesX}

\description{
  Run BayesX program files from \R.
}

\usage{
run.bayesx(dir, prg.name = "bayesx.estim.input.prg", 
  verbose = FALSE, bin = getOption("bayesx.bin"))
}

\arguments{
  \item{dir}{a character string, the directory where the prgram file is stored.}
  \item{prg.name}{the name of the program file.}
  \item{verbose}{should output be printed to the \R console during runtime of 
    \code{\link{bayesx}}.}
  \item{bin}{set the \code{bin} path where the executable binary of \pkg{BayesX} command-line 
    version is installed, may also be set in function \code{\link{options}} via 
    \code{options(bayesx.bin="/path/to/BayesX")}, where \code{"/path/to/BayesX"} is the bin path 
    provided as character string.}
}

\details{
  Function uses \code{\link[base]{system}} to run BayesX within an \R session.
}

\value{
  Returns a \code{list} containg information if \pkg{BayesX} was succesfully launched and how long 
  the process was running.
}

\author{
  Thomas Kneib, Stefan Lang, Nikolaus Umlauf, Achim Zeileis.
}

\section{WARNINGS}{
  Function \code{\link{run.bayesx}} cannot be processed if the \code{bin} path to the BayesX binary 
  executable command-line version is not specified properly, i.e. before using \code{run.bayesx} the 
  binary has to be installed, see function \code{\link{install.bayesx}} for installation details. 
  After installation the \code{bin} path can bet set with

  \code{options(bayesx.bin = "/path/to/BayesX")}

  where \code{"/path/to/BayesX"} represents the full path to the binary with the name of the binary 
  at last, e.g. on Windows systems this may be something like 
  \code{"C:/programs/BayesX/commandline/bayesx.exe"}. 
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

## first create BayesX .prg file
pars <- parse.bayesx.input(y ~ s(x, bs = "ps") + 
  s(z, w, bs = "te") + fac, data = dat)
prg <- write.bayesx.input(pars)

## now estimate model with BayesX
ok <- run.bayesx(prg$file.dir)
writeLines(ok$log)
print(ok$runtime)
}
}

\keyword{regression}
