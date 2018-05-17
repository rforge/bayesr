\name{read.bnd}
\alias{read.bnd}

\title{Read Geographical Information in Boundary Format}

\description{
Reads the geographical information provided in a file in boundary format (see Ch. 5 of the BayesX Reference Manual) and
stores it in a map object.
}

\usage{
read.bnd(file, sorted=FALSE)
}

\arguments{
  \item{file}{Name of the boundary file to be read.}
  \item{sorted}{Should the regions be ordered by the numbers speciying the region names (\code{sorted=TRUE})?}
}

\value{Returns a list of polygons that form the map. Additional attributes are
  \item{surrounding}{Parallel list where for each polygon, the name of a possible
  surrounding region is saved.}
  \item{height2width}{Ratio between height and width of the map. Allows customised drawing and storage in files by specifying the appropriate height and width.}
  \item{class}{Indicates whether the map is stored in boundary format (\code{bnd}) or graph format (\code{gra}). Maps returned by \code{read.bnd} are of class \code{bnd}}
}

\references{BayesX Reference Manual. Available from \url{http://www.bayesx.org/}}

\author{Daniel Sabanes Bove, Felix Heinzl, Thomas Kneib, Andreas Brezger}

\seealso{\code{\link{write.bnd}},\code{\link{drawmap}},\code{\link{read.gra}},\code{\link{write.gra}}}

\examples{
germany <- read.bnd(system.file("examples/germany.bnd", package="BayesX"))
drawmap(map=germany)
attributes(germany)

germany <- read.bnd(system.file("examples/germany2001.bnd", package="BayesX"))
drawmap(map=germany)
attributes(germany)
}

\keyword{spatial}
