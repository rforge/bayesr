\name{smooth.construct}
\alias{smooth.construct}

\title{Constructor Functions for Smooth Terms in BAMLSS}

\description{The generic function is only a copy of \code{\link[mgcv]{smooth.construct}} adding
  a \code{\dots} argument. For objects of class \code{"bamlss.frame"} and \code{"bamlss"} the
  method extracts all smooth model terms, see function \code{\link{bamlss.frame}} for details
  on the setup of BAMLSS.}

\usage{
smooth.construct(object, data, knots, ...)
}

\arguments{
  \item{object}{Either a smooth specifiaction object, or object of class \code{"bamlss"} or
    \code{"bamlss.frame"}. For smooth specification objects, see function
    \code{\link[mgcv]{smooth.construct}}}
  \item{data}{A data frame or list, see also see function \code{\link[mgcv]{smooth.construct}}}
  \item{knots}{See function \code{\link[mgcv]{smooth.construct}}}
  \item{\dots}{Arguments passed to the smooth term constructor functions.}
}

\value{
For smooth specification objects see function see \code{\link[mgcv]{smooth.construct}}. For
objects of class \code{"bamlss.frame"} or \code{"bamlss"} the list of smooth constructs,
see function \code{\link{bamlss.frame}} for more details.
}

\seealso{
  \code{\link{bamlss.frame}}, \code{\link{bamlss}}, \code{\link[mgcv]{smooth.construct}}.
}


\examples{
## Generate some data.
d <- GAMart()

## Create a "bamlss.frame".
bf <- bamlss.frame(num ~ s(x1) + s(x2), data = d)

## Extract the smooth construct.
sc <- smooth.construct(bf)
str(sc)

## Also possible with formulas.
f <- bamlss.formula(list(
  num ~ s(x1) + te(lon,lat),
  sigma ~ s(x2)
), family = "gaussian")

sc <- smooth.construct(f, data = d)
str(sc)
}

\keyword{models} \keyword{smooth} \keyword{regression}

