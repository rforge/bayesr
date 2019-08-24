\name{fuse}
\alias{fuse}

\title{Combine Regions}

\description{
Combines a list of several regions of a map object in boundary format into a single region.
}

\usage{
fuse(map, regions, name)
}

\arguments{
  \item{map}{Map object in boundary format that should be modified.}
  \item{regions}{Vector of regions to be combined}
  \item{name}{Name that should be given to the region arising from fusing the specified regions.}
}

\value{
Map object in boundary format with the specified regions combined.
}

\author{Nadja Klein}

\seealso{\code{\link{read.bnd}},\code{\link{write.bnd}}}

\examples{
\dontrun{require("gpclib")
library("maptools")
gpclibPermit()
map <- read.bnd(system.file("examples/germany9301.bnd", package="BayesX"))
drawmap(map=map, drawnames=TRUE)

#vector of regions to be combined 
regions <- c("1056","1060","1061")
#new name of combined region
newname <- "1"
newmap <- fuse(map,regions,newname)
drawmap(map=newmap,drawnames=TRUE)

#vector of regions to be combined
germany <- read.bnd(system.file("examples/germany.bnd", package="BayesX"))
drawmap(map=germany, drawnames=TRUE)
regions <- c("9371","9373","9374","9471","9472","9474","9574")
#new name of combined region
newname <- "1"
newmap <- fuse(germany,regions,newname)
drawmap(map=newmap,drawnames=TRUE)}
}


\keyword{spatial}
