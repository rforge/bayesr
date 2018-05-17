\name{shp2bnd}
\alias{shp2bnd}

\title{convert a shape-file into a boundary object}

\description{
Converts the geographical information provided in a shape-file into a boundary object
(see Ch. 5 of the Reference Manual)
}

\usage{
shp2bnd(shpname, regionnames, check.is.in = TRUE)
}

\arguments{
  \item{shpname}{Base filename of the shape-file (including path)}
  \item{regionnames}{Either a vector of region names or the name of the variable in the dbf-file representing these names}
  \item{check.is.in}{Test whether some regions are surrounded by other regions (FALSE speeds up the execution time but may result in a corrupted bnd-file)}
}

\value{Returns a boundary object, i.e. a list of polygons that form the
  map. See \code{\link{read.bnd}} for more information on the format.
}

\references{BayesX Reference Manual. Available from \url{http://www.bayesx.org}}

\author{Felix Heinzl, Daniel Sabanes Bove, Thomas Kneib with contributions by Michael Hoehle and Frank Sagerer}

\seealso{\code{\link{write.bnd}},\code{\link{drawmap}},\code{\link{read.bnd}}}

\examples{
## read shapefile into bnd object
shpName <- sub(pattern="(.*)\\\\.dbf", replacement="\\\\1",
               x=system.file("examples/northamerica_adm0.dbf",
                             package="BayesX")) 
north <- shp2bnd(shpname=shpName, regionnames="COUNTRY")

## draw the map
drawmap(map=north)

## compare with shipped bnd file
shippedBnd <- read.bnd(system.file("examples/northamerica.bnd", package="BayesX"))
stopifnot(all.equal(north, shippedBnd))
}

\keyword{spatial}
