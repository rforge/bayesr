\name{bnd2gra}
\alias{bnd2gra}

\title{Convert Boundary Format to Graph Format}

\description{
Converts a map in boundary format to a map in graph format.
}

\usage{
bnd2gra(map)
}

\arguments{
  \item{map}{Map in boundary format that should be converted.}
}

\value{
Returns an adjacency matrix that represents the neighborhood structure of the map object in graph format.
}

\references{BayesX Reference Manual. Available from \url{http://www.bayesx.org}}

\author{Felix Heinzl, Thomas Kneib}

\seealso{\code{\link{read.bnd}},\code{\link{read.gra}},\code{\link{write.bnd}},\code{\link{write.gra}}}

\examples{
tanzania.bnd <- read.bnd(system.file("examples/tanzania.bnd", package="BayesX"))
tanzania.gra <- bnd2gra(tanzania.bnd)
}

\keyword{spatial}
