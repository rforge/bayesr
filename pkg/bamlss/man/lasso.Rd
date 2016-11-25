\name{la}
\alias{la}
\alias{lasso}

\title{Lasso Smooth Constructor}

\description{
  Smooth constructors for Lasso penalization with \code{\link{bamlss}}. The penalization is based
  on a Taylor series approximation of the Lasso penalty.
}

\usage{
la(formula, type = c("single", "multiple"), ...)
}

\arguments{
  \item{formula}{A formula like \code{~ x1 + x2 + ... + xk} of variables which should be
    penalized with Lasso.}
  \item{type}{Should one single penalty parameter be used or multiple parameters, one for each
    covariate in \code{formula}.}
  \item{\dots}{Arguments passed to the subsequent smooth constructor function.
    \code{lambda} controls the starting value of the penalty parameter, \code{const} the constant
    that is added within the penalty approximation.}
}

\value{
  Similar to function \code{\link[mgcv]{s}} a simple smooth specification object.
}

\references{
  Oelker Margreth-Ruth and Tutz Gerhard (2015). A uniform framework for combination of
  penalties in generalized structured models. \emph{Adv Data Anal Classif}.
  \url{http://dx.doi.org/10.1007/s11634-015-0205-y}
}

\seealso{
\code{\link[mgcv]{s}}, \code{\link[mgcv]{smooth.construct}}
}

\examples{
\dontrun{## Simulate data.
dgp <- function(n = 100, p = 10,
  sd = 0.3, sigma = diag(p), prob = c(0.8, 0.2), kmax = NULL)
{
  require("mvtnorm")

  beta <- if(is.null(kmax)) {
    sample(0:1, size = p, replace = TRUE, prob = prob)
  } else {
    sample(c(rep(0, p - kmax), rep(1, kmax)))
  }

  X <- rmvnorm(n, rep(0, p), sigma = sigma)
  colnames(X) <- paste("x", 1:ncol(X), sep = "")
  f <- X %*% beta
  y <- f + rnorm(n, sd = sd)

  d <- data.frame("y" = as.vector(y), "f" = as.vector(f))
  d <- cbind(d, as.data.frame(X))

  attr(d, "beta") <- beta
  return(d)
}

set.seed(123)
d <- dgp(n = 1000, p = 20, sd = 0.1)

f <- grep("x", colnames(d), value = TRUE)
f <- paste(f, collapse = "+")

f0 <- as.formula(paste("y", f, sep = "~"))
f <- as.formula(paste("", f, sep = "~"))

b0 <- bamlss(f0, data = d, sampler = FALSE)
b1 <- bamlss(y ~ la(f), data = d, sampler = FALSE)
b2 <- bamlss(y ~ la(f, type = "multiple"), data = d, sampler = FALSE)

cb0 <- coef(b0, hyper = FALSE)
cb0 <- cb0[!grepl("Intercept", names(cb0))]

cb1 <- coef(b1, hyper = FALSE)
cb1 <- cb1[!grepl("Intercept", names(cb1))]

cb2 <- coef(b2, hyper = FALSE)
cb2 <- cb2[!grepl("Intercept", names(cb2))]

coefs <- data.frame(
  "true" = attr(d, "beta"),
  "fixed" = cb0,
  "lasso" = cb1,
  "sep.lasso" = cb2
)

print(round(coefs, 5))
}
}

\keyword{regression}

