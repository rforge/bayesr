\name{read.gra}
\alias{read.gra}

\title{Read Geographical Information in Graph Format}

\description{
Reads the geographical information provided in a file in graph format (see Ch. 5 of the BayesX Reference Manual) and
stores it in a map object.
}

\usage{
read.gra(file, sorted=FALSE)
}

\arguments{
  \item{file}{Name of the graph file to be read.}
  \item{sorted}{Should the regions be ordered by the numbers speciying the region names (\code{sorted=TRUE})?}
}

\value{
Returns an adjacency matrix that represents the neighborhood structure defined in the graph file. Additional attributes are
  \item{dim}{Dimension of the (square) adjacency matrix.}
  \item{dimnames}{List of region names corresponding to rows and columns of the adjacency matrix.}
  \item{class}{Indicates whether the map is stored in boundary format (\code{bnd}) or graph format (\code{gra}). Maps returned by \code{read.gra} are of class \code{gra}}
}

\references{BayesX Reference Manual. Available from \url{http://www.bayesx.org/}}

\author{Thomas Kneib, Felix Heinzl}

\seealso{\code{\link{write.gra}},\code{\link{read.bnd}},\code{\link{write.bnd}},\code{\link{get.neighbor}},\code{\link{add.neighbor}},\code{\link{delete.neighbor}}}

\examples{
germany <- read.gra(system.file("examples/germany.gra", package="BayesX"))
attributes(germany)
}

\keyword{spatial}
