\name{install.bayesx}
\alias{install.bayesx}

\title{
  Install BayesX
}

\description{
  Function to Install \pkg{BayesX}.
}

\usage{
install.bayesx(inst.dir = NULL, source.dir = NULL)
}

\arguments{
  \item{inst.dir}{a character string, specifies the directory \pkg{BayesX} should be installed in.}
  \item{source.dir}{character string, the directory where on windows systems the installation file 
    \code{BayesX_windows.exe} lies in, else the file \code{bayesxsource.zip}.}
}

\details{
  On windows platforms function \code{install.bayesx} tries to install the file 
  \code{BayesX_windows.exe}, on unix platforms the \pkg{BayesX} sources in file 
  \code{bayesxsource.zip} are tried to compiled with the GNU C++ compiler. On all platforms, if 
  \code{source.dir = NULL} it is tried to download the necessary installation files from 
  \url{http://www.stat.uni-muenchen.de/~bayesx/bayesx.html}, and if \code{inst.dir = NULL}, 
  function \code{install.bayesx} tries to install \pkg{BayesX} in the \code{\link{.libPaths}} 
  directory this package lies in.
}

\note{
  During the installation process on Windows systems, you will be asked to locate java.exe for 
  BayesX.
}

\value{
  If \pkg{BayesX} is successfully installed, a character string of the installation path of the 
  command-line version executable binary is returned.
}

\author{
  Thomas Kneib, Stefan Lang, Nikolaus Umlauf, Achim Zeileis.
}

\seealso{
  \code{\link{check.install.bayesx}}.
}

\examples{
\dontrun{
install.bayesx()
}
}

\keyword{regression}
