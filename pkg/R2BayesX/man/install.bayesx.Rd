\name{install.bayesx}
\alias{install.bayesx}

\title{
  Install BayesX
}

\description{
  Function to Install the BayesX command-line binary.
}

\usage{
install.bayesx(inst.dir = NULL, source.dir = NULL, type = NULL)
}

\arguments{
  \item{inst.dir}{a character string, specifies the directory \pkg{BayesX} should be installed in.}
  \item{source.dir}{character string, the directory where on Windows systems the installation file 
    \code{BayesX_windows.exe} lies in, else the file \code{bayesxsource.zip}.}
  \item{type}{to force downloading the BayesX sources \code{bayesxsource.zip} and compile with GNU
    Compiler Collection (GCC) C compiler, set \code{type = "source"}.}
}

\details{
  On windows platforms function \code{install.bayesx} tries to install the file 
  \code{BayesX_windows.exe}, on unix platforms the BayesX sources in file \code{bayesxsource.zip} 
  are tried to compiled with the GNU C++ compiler. On all platforms, if \code{source.dir = NULL} it 
  is tried to download the necessary installation files from 
  \url{http://www.stat.uni-muenchen.de/~bayesx/bayesx.html}, and if \code{inst.dir = NULL}, function 
  \code{install.bayesx} tries to install the BayesX binary in the \code{\link{.libPaths}} 
  directory this package is saved.
}

\note{
  During the installation process on Windows systems, you will be asked to locate java.exe for 
  BayesX.
}

\value{
  If the BayesX binary is successfully installed, a character string of the installation path of the 
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
