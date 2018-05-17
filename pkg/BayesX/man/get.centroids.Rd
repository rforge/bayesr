\name{get.centroids}
\alias{get.centroids}

\title{Compute Centroids of Polygons}

\description{
Computes all areas and centroids of the regions of a given map in boundary format.
}

\usage{
get.centroids(map)
}

\arguments{
  \item{map}{Map object in boundary format.}
}

\value{
Matrix of area and centroids.
}

\author{Felix Heinzl, Thomas Kneib}

\examples{
germany <- read.bnd(system.file("examples/germany.bnd", package="BayesX"))
centroids <- get.centroids(germany)
centroids[1:10,]

plot(c(2100,3700),c(6800,8500),type="n", xlab="", ylab="")
for(i in 1:10){
   polygon(germany[[i]])
   region <- attr(germany,"names")[i]
   text(x=centroids[i,2]+50, y=centroids[i,3]+30, region, cex=0.7)
}   
points(centroids[1:10,2:3], col='red', pch=16)
}

\keyword{spatial}
