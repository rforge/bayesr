\name{write.gra}
\alias{write.gra}

\title{Saving Maps in Graph Format}

\description{
Writes the information of a map object to a file (in graph format)
}

\usage{
write.gra(map, file, replace = FALSE)
}

\arguments{
  \item{map}{Map object ot be saved (should be in graph format, see \code{\link{bnd2gra}} for the conversion of boundary format to graph format).}
  \item{file}{Name of the file to write to}
  \item{replace}{Should an existing file be overwritten with the new version?}
}

\references{BayesX Reference Manual. Available from \url{http://www.bayesx.org}}

\author{Thomas Kneib, Felix Heinzl}

\seealso{\code{\link{write.bnd}},\code{\link{read.gra}},\code{\link{read.bnd}}}

\keyword{spatial}
