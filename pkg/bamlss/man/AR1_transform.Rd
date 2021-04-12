\name{AR1_transform}
\alias{AR1_transform}

\title{
  AR1 Transformation
}

\description{
  The transformer function takes a \code{\link{bamlss.frame}} object and transforms
  the response and the design matrices to account for lag 1 autocorrelation.
}

\usage{
AR1_transform(x)
}

\arguments{
  \item{x}{Object returned from function \code{\link{bamlss.frame}}.}
}

\value{
  A transformed \code{\link{bamlss.frame}}.
}

\seealso{
  \code{\link{bamlss.frame}}, \code{\link{bamlss}}, \code{\link[mgcv]{smooth2random}}.
}

\examples{\dontrun{## Simulate data.
set.seed(123)

n <- 1000
d <- data.frame("x" = runif(n, -3, 3))
d$y <- rep(0, n)
for(i in 2:n) {
  d$y[i] <- sin(d$x[i]) + 0.3 * d$y[i - 1] + rnorm(1, sd = 0.1)
}

plot(d)

b <- bamlss(y ~ s(x), data = d, transform = AR1_transform(rho = 0.3))

e <- residuals(b, type = "response")
acf(e)
}}

\keyword{regression}
\keyword{smooth}

