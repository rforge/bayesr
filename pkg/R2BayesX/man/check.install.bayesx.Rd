\name{check.install.bayesx}
\alias{check.install.bayesx}

\title{
  Check if Bayesx is Installed
}

\description{
  This function simply checks if BayesX command-line binary is already installed and ready for use.
}

\usage{
check.install.bayesx(bin = getOption("bayesx.bin"), verbose = FALSE)
}

\arguments{
  \item{bin}{the full path to the executable command-line version binary of \pkg{BayesX}, may also 
    be set in function \code{\link{options}} via \code{options(bayesx.bin = "/path/to/BayesX")}, where 
    \code{"/path/to/BayesX"} is the \code{bin} path provided as character string, e.g. on windows 
   systems this may be something like \code{"C:/programs/BayesX/commandline/bayesx.exe"}.}
  \item{verbose}{if set to \code{TRUE} some additional information will be printed on the console.}
}

\details{
  Function \code{check.install.bayesx} tries to execute a mini program in \pkg{BayesX} to verify the 
  installation of the command-line version executable binary.
}

\value{
  If the BayesX binary is successfully installed, a character string of the full installation path 
  of the command-line version executable binary is returned, if the path to BayesX is already set in 
  the environment PATH variable of the operating system, only the name of the binary is returned, in 
  all other cases \code{NULL} is returned.
}

\author{
  Thomas Kneib, Stefan Lang, Nikolaus Umlauf, Achim Zeileis.
}

\seealso{
  \code{\link{install.bayesx}}.
}

\examples{
check.install.bayesx()
}

\keyword{regression}
