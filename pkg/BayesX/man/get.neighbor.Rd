\name{get.neighbor}
\alias{get.neighbor}

\title{Obtain Neighbors of Given Regions}

\description{
Extracts the neighbors of a number of regions from a map in graph format.
}

\usage{
get.neighbor(map, regions)
}


\arguments{
  \item{map}{Map object in graph format.}
  \item{regions}{Vector of names of regions for which the neighbors should be axtracted.}
}

\value{
A list of vectors containing the neighbors of the elements in \code{regions}.
}

\author{Felix Heinzl, Thomas Kneib}

\seealso{\code{\link{add.neighbor}},\code{\link{delete.neighbor}}}

\examples{
germany <- read.gra(system.file("examples/germany.gra", package="BayesX"))
get.neighbor(germany, "1001")
get.neighbor(germany, c("1001", "7339"))
}

\keyword{spatial}
