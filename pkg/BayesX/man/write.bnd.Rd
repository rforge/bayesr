\name{write.bnd}
\alias{write.bnd}

\title{Saving Maps in Boundary Format}

\description{
Writes the information of a map object to a file (in boundary format)
}

\usage{
write.bnd(map, file, replace = FALSE)
}

\arguments{
  \item{map}{Map object ot be saved (should be in boundary format).}
  \item{file}{Name of the file to write to}
  \item{replace}{Should an existing file be overwritten with the new version?}
}

\references{BayesX Reference Manual. Available from \url{http://www.bayesx.org}}

\author{Thomas Kneib, Felix Heinzl}

\seealso{\code{\link{write.gra}},\code{\link{read.gra}},\code{\link{read.bnd}}}

\keyword{spatial}
