\name{plotsurf}
\alias{plotsurf}
\alias{akimaPermit}
\alias{akimaPermitStatus}

\title{Visualise Surface Estimates}
\description{
Visualises surface estimates obtained with BayesX.
}
\usage{
plotsurf(data, x=2, y=3, z=4, mode=1, ticktype="detailed", 
         expand=0.75, d=100, theta=-30, phi=25, ...)
akimaPermitStatus()
akimaPermit()
}

\arguments{
  \item{data}{Either the name of a file or a data frame containing the estimation results.}
  \item{x}{Defines the x-axis in the plot. Either the name of a variable in \code{data} or the index of the corresponding column.}
  \item{y}{Defines the y-axis in the plot. Either the name of a variable in \code{data} or the index of the corresponding column.}
  \item{z}{Defines the z-axis in the plot. Either the name of a variable in \code{data} or the index of the corresponding column.}
  \item{mode}{\code{plotsurf} is mostly a convenient interface to the functions \code{persp} (\code{mode=1}), \code{image} (\code{mode=2}) and
        \code{contour} (\code{mode=3}).}
  \item{ticktype, expand, d, theta, phi}{Overwrite the default behaviour of \code{persp}}
  \item{\dots}{Further parameteres that are parsed to the internal call to \code{persp}, \code{image} or \code{contour}}
}

\details{The akimaPermit function is used to choose to permit the use of akima if installed,
and akimaPermitStatus reports its status. The licence for akima is not Free or Open Source and explicitly forbids commercial use.
See the akima licence file for details.}

\author{Felix Heinzl, Thomas Kneib}

\seealso{\code{\link{drawmap}},\code{\link{plotautocor}},\code{\link{plotsample}},\code{\link{plotnonp}}}

\examples{
res <- read.table(system.file("examples/surface_f_x1_x2_pspline.res", 
                              package="BayesX"), header=TRUE)

# call akimaPermit to allow use of akima despite its licence that restricts
# usage to non-commercial purposes
akimaPermit()
plotsurf(res)
plotsurf(res, mode=2)
plotsurf(res, mode=3)

plotsurf(res, x="x1", y="x2", z="pmed")

plotsurf(res, ticktype="simple")

plotsurf(res, main="3D-Plot", xlab="myx", ylab="myy", zlab="f(myx,myy)")
}

\keyword{hplot}
