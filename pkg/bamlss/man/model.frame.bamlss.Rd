\name{model.frame.bamlss}
\alias{model.frame.bamlss}
\alias{model.frame.bamlss.frame}

\title{
  BAMLSS Model Frame
}

\description{
  Extracts the model frame of a \code{\link{bamlss}} or \code{\link{bamlss.frame}}
  object.
}

\usage{
\method{model.frame}{bamlss}(formula, ...)

\method{model.frame}{bamlss.frame}(formula, ...)
}

\arguments{
  \item{formula}{An object of class \code{"bamlss"} or \code{"bamlss.frame"}.}
  \item{\dots}{Other arguments.}
}

\value{
  The \code{\link[base]{data.frame}} containing the variables used for modeling
  with function \code{\link{bamlss}}.
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}},
  \code{\link{model.matrix.bamlss.frame}}.
}

\examples{
\dontrun{## Generate some data.
d <- GAMart()

## Model formula.
f <- list(
  log(pnum) ~ s(x1) + log(x2),
  sigma ~ s(x3)
)

## Estimate model.
b <- bamlss(f, data = d)

## Extract the model frame.
head(model.frame(b))
}
}

\keyword{regression}
\keyword{models
}
