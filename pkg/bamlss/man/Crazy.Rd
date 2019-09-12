\name{Crazy}
\alias{Crazy}

\title{
  Crazy simulated data
}

\description{
  This function creates simulated data based on a crazy function.
}

\usage{
Crazy(n)
}

\arguments{
  \item{n}{The number of observations to be simulated.}
}

\value{
  A data frame with response \code{y} and covariate \code{x}.
}

\seealso{
\code{\link{GAMart}}
}

\examples{
d <- Crazy()
head(d)
plot(d)
}

\keyword{datasets}

