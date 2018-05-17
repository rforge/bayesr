\name{smooth.bnd}
\alias{smooth.bnd}

\title{Round Boundary Information}

\description{
Rounds the boundary information in a map object in boundary format to a specified precision.
}

\usage{
smooth.bnd(map, digits = 2, scale = 1)
}

\arguments{
  \item{map}{Map object in boundary format that should be modified.}
  \item{digits}{Number of digits to round to.}
  \item{scale}{Scaling factor that should be applied for rounding. For example, with
               \code{scale=0.1} all polygons are divided by 10 before rounding.}
}

\value{
Map object in boundary format rounded to the specified precision.
}

\author{Felix Heinzl, Thomas Kneib}

\seealso{\code{\link{read.bnd}},\code{\link{write.bnd}}}

\keyword{spatial}
