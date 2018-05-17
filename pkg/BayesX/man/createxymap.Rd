\name{createxymap}
\alias{createxymap}

\title{Create map objects for some points and a given distance}

\description{
Creates a map object from a list of coordinates by treating observations within a certain distance as neighbors. 
The resulting neighborhood structure is stored in a map object in graph format while a map
in boundary format is created to enable visualisation.
}

\usage{
createxymap(x, y, districts=NULL, p=2, max.dist)
}

\arguments{
  \item{x}{Vector of x-coordinates of underlying points}
  \item{y}{Vector of y-coordinates of underlying points}
  \item{districts}{Either \code{NULL} or a vector of names for labeling points. If \code{districts=NULL}, points are labelled by index.}
  \item{p}{Any p-norm with p>=1 can be chosen as the distance measure with the Euclidian distance (\code{p=2}) being the default. For \code{p=Inf}, the maximum of coordinates is used. Except for \code{p=Inf}, calculations can fail for huge \code{p}.}
  \item{max.dist}{Value which determines the neighborship. Points whose distance is smaller or equal than \code{max.dist} are considered as neighbors.}
}

\author{Felix Heinzl, Thomas Kneib}

\value{
List of two elements: map object in graph format and map object in boundary format.
}

\seealso{\code{\link{read.gra}},\code{\link{read.bnd}},\code{\link{drawmap}}}

\examples{
x <- c(3,3,5,9.5,11,11)
y <- c(2,4,4,6,4.5,5)
xymap <- createxymap(x, y, districts=c("A","B","C","D","E","F"), max.dist=2)
xymap$gra
drawmap(map=xymap$bnd)
}

\keyword{spatial}
