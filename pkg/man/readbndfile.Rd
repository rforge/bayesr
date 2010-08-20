\name{readbndfile}
\alias{readbndfile}

\title{Read a Map Boundary File}
\description{This functions reads in a \code{\link{polygon}} boundary file of a map}

\usage{
readbndfile(path = path, name = name)
}

\arguments{
  \item{path}{The path of the boundary file}
  \item{name}{The name of the new generated map object}
}
 
\value{A map object with a \code{\link{list}} of \code{\link{polygon}}s is returned and stored in the object specified
       by argument \code{name}}

\author{Thomas Kneib \email{thomas.kneib@stat.uni-muenchen.de}}

\note{For further details on map and shapefile objects see \url{http://www.stat.uni-muenchen.de/~bayesx/bayesxmaps.html}, or
      \url{http://www.stat.uni-muenchen.de/~kneib/software.html}}

\seealso{\code{\link{drawmap}}}

\examples{
# read in a boundary file of munich districts
readbndfile("http://www.stat.uni-muenchen.de/~bayesx/maps/munich_bezt.bnd","munich")
drawmap(munich,main="District map of Munich, Germany")

# map of nigeria
readbndfile("http://www.stat.uni-muenchen.de/~bayesx/maps/nigeria.bnd","nigeria")
drawmap(nigeria,main="Map of Nigeria")

# map of germany
readbndfile("http://www.stat.uni-muenchen.de/~bayesx/maps/germany2001.bnd","germany")
drawmap(germany,main="Map of Germany")
}

\keyword{aplot}
