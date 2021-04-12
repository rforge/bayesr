\name{AR1_transform}
\alias{AR1_transform}

\title{
  AR1 Transformer Function
}

\description{
  The transformer function takes a \code{\link{bamlss.frame}} object and transforms
  the response and the design matrices to account for lag 1 autocorrelation. The method
  is also known as Prais-Winsten estimation.
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

\references{
  Johnston, John (1972). Econometric Methods (2nd ed.). New York: McGraw-Hill. pp. 259--265.
}

\seealso{
  \code{\link{bamlss.frame}}, \code{\link{bamlss}}, \code{\link[mgcv]{smooth2random}}.
}

\examples{\dontrun{## Simulate AR1 data.
set.seed(123)

n <- 1000
d <- data.frame("t" = 1:n)
f <- function(x) {
  sin(x / n * 2 * pi - pi)
}
d$y <- rep(0, n)
for(i in 2:n) {
  d$y[i] <- f(d$t[i]) + 0.9 * (d$y[i - 1] - f(d$t[i - 1])) + rnorm(1, sd = 0.05)
}

plot(d, type = "l")

## Estimate models without and with AR1 transformation.
b0 <- bamlss(y ~ s(t,k=40), data = d)
b1 <- bamlss(y ~ s(t,k=40), data = d, transform = AR1_transform(rho = 0.9))

## Estimate AR1 correlation.
b2 <- bamlss(y ~ s(t,k=40), data = d, family = "AR1")
rho <- predict(b2, model = "alpha", type = "parameter")
print(rho[1])

## Plot fitted trends.
p0 <- predict(b0, model = "mu")
p1 <- predict(b1, model = "mu")
p2 <- predict(b2, model = "mu")

plot(d, type = "l")
lines(f(d$t) ~ d$t, col = 2, lwd = 2)
lines(p0 ~ d$t, col = 4, lwd = 2)
lines(p1 ~ d$t, col = 3, lwd = 3)
lines(p2 ~ d$t, col = 5, lwd = 3)
legend("topleft", c("no trans", "with trans", "AR model", "truth"),
  lwd = 2, col = c(4, 3, 5, 2), bty = "n")
}}

\keyword{regression}
\keyword{smooth}

