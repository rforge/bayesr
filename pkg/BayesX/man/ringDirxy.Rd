\name{.ringDirxy}
\alias{.ringDirxy}

\title{Compute direction of points in a polygon.}
\description{
  This is (private, not exported) R code replacing the call to
  maptools:::.ringDirxy, which links to C code.
}
\usage{
.ringDirxy(xy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xy}{a two-column matrix with x and y coordinates of at least 3 points which form a polygon}
}
\value{
  Returns the integer -1 if the direction is counter-clockwise and 1 if
  the direction is clockwise.
}
\note{
  The function uses vectorization and should not be very much slower
  than the original.
}
\author{Daniel Sabanes Bove}

% we don't run the functions to not depend on the private
% function in maptools
\examples{
## create test polygon ordered clockwise
nPoints <- 200
sequence <- seq(from=1, to=10, length=nPoints)
firstHalf <- sequence[1:(nPoints/2)]

test <- cbind(x=c(sequence, rev(sequence)),
              y=c(firstHalf, rev(firstHalf), -firstHalf, -rev(firstHalf)))
plot(test,
     type="n")
polygon(test)

## compute direction with this and the original function
stopifnot(identical(BayesX:::.ringDirxy(test), as.integer(+1)))
\dontrun{maptools:::.ringDirxy(test)}

## now transform to counter-clockwise and try that
test[, 2] <- - test[, 2]
stopifnot(identical(BayesX:::.ringDirxy(test), as.integer(-1)))
\dontrun{maptools:::.ringDirxy(test)}
}

\keyword{internal}

