\name{family.bamlss}
\alias{family.bamlss}
\alias{family.bamlss.frame}
\alias{bamlss.family}
\alias{beta.bamlss}
\alias{binomial.bamlss}
\alias{cnorm.bamlss}
\alias{cox.bamlss}
\alias{gaussian.bamlss}
\alias{gamma.bamlss}
\alias{multinomial.bamlss}
\alias{mvnorm.bamlss}
\alias{poisson.bamlss}
\alias{jm.bamlss}

\title{
  Distribution Families in \pkg{bamlss}
}

\description{
  Family objects in \pkg{bamlss} specify the information that is needed for using (different)
  model fitting engines, e.g., the parameter names and corresponding link functions, the
  density function, derivatives of the log-likelihood w.r.t. the linear predictors, and so
  forth. The optimizer or sampler functions that are called by \code{\link{bamlss}} must know how
  much information is needed to interpret the model since the family objects are simply
  passed through. Family objects are also used for computing post-modeling statistics, e.g., for
  residual diagnostics or random number generation. See the details and examples.
}

\usage{
## Family objects in bamlss:
beta.bamlss(...)
binomial.bamlss(...)
cnorm.bamlss(...)
cox.bamlss(...)
jm.bamlss(...)
gaussian.bamlss(...)
gamma.bamlss(...)
multinomial.bamlss(...)
mvnorm.bamlss(...)
poisson.bamlss(...)

## Extractor functions:
\method{family}{bamlss}(object, ...)
\method{family}{bamlss.frame}(object, ...)
}

\arguments{
  \item{object}{An object of class \code{"bamlss"} or \code{"bamlss.frame"}, see function
    \code{\link{bamlss}} and \code{\link{bamlss.frame}}.}
  \item{\dots}{Arguments passed to functions that are called within the family object.}
}

\details{
  The following lists the minimum requirements on a \pkg{bamlss} family object to be used with
  \code{\link{bamlss}} and \code{\link{bamlss.frame}}:
  \itemize{
    \item The family object must return a \code{\link{list}} of class \code{"family.bamlss"}.
    \item The object must contain the family name as a character string. 
    \item The object must contain the names of the parameters as a character string, as well as
      the corresponding link functions as character string.
  }
  For most optimizer and sampling functions at least the density function should be provided.
  When using generic model fitting engines like \code{\link{bfit}} or \code{\link{GMCMC}}, as
  well as for computing post-modeling statistics with function \code{\link{samplestats}}, and others,
  it is assumed that the density function in a family object has the following arguments:

  \code{d(y, par, log = FALSE, ...)}

  where argument \code{y} is the response (possibly a matrix) and \code{par} is a named list
  holding the evaluated parameters of the distribution, e.g., using a normal distribution \code{par}
  has two elements, one for the mean \code{par$mu} and one for the standard deviation
  \code{par$sigma}. The dots argument is for passing special internally used objects, depending
  on the type of model this feature is usually not needed.

  Similarly, for derivative based algorithms, e.g. using iteratively reweighted least squares
  (IWLS, see function \code{\link{bfit}}, the family object holds derivative functions evaluating
  derivatives of the log-likelihood w.r.t. the linear predictors (or expectations of derivatives).
  For each parameter, these functions also hold the following arguments:
 
  \code{score(y, par, ...)}

  for computing the first derivative of the log-likelihood w.r.t. a linear predictor and

  \code{hess(y, par, ...)}

  for computing the negative second derivatives. Within the family object these functions
  are organized in a named list, see the examples below.

  In addition, for the cumulative distribution function (\code{p(y, par, ...)}), for the quantile
  function (\code{q(y, par, ...)}) or for creating random numbers (\code{r(y, par, ...)}) the same
  structure is assumed. See, e.g., the code of function \code{gaussian.bamlss()}.

  Some model fitting engines can initialize the distributional parameters which oftentimes
  leads to much faster convergence. The initialize functions are again organized witin a named list,
  one entry for each parameter, similar to the \code{score} and \code{hess} functions, e.g.,
  see the code of family object \code{gaussian.bamlss}.

  Using function \code{\link{bamlss}}, \code{\link{residuals.bamlss}} and
  \code{\link{predict.bamlss}} the family objects may
  also specify the \code{transform()}er, \code{optimizer()}, \code{sampler()},
  \code{samplestats()}, \code{results()}, \code{residuals()} and \code{predict()} function that
  should be used with this family. See for example the setup of \code{\link{cox.bamlss}}.

  For using specialized estimation engines like \code{\link{JAGS}} it is recommended to supply
  any extra arguments needed by those engines with an additional list entry within the family
  object, e.g., using \code{gaussian.bamlss} with \code{\link{JAGS}} the family objects holds
  special details in an element named \code{"bugs"}.

  The examples below illustrate this setup. See also the code of the \pkg{bamlss} family
  functions.
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}}
}

\examples{
## New family object for the normal distribution,
## can be used by function bfit() and GMCMC().
normal.bamlss <- function(...) {
  f <- list(
    "family" = "normal",
    "names" = c("mu", "sigma"),
    "links" = c("identity", "log"),
    "d" = function(y, par, log = FALSE) {
      dnorm(y, mean = par$mu, sd = par$sigma, log = log)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        drop((y - par$mu) / (par$sigma^2))
      },
      "sigma" = function(y, par, ...) {
        drop(-1 + (y - par$mu)^2 / (par$sigma^2))
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        drop(1 / (par$sigma^2))
      },
      "sigma" = function(y, par, ...) { 
        rep(2, length(y))
      }
    )
  )
  class(f) <- "family.bamlss"
  return(f)
}

\dontrun{## Test on simulated data.
d <- GAMart()
b <- bamlss(num ~ s(x1) + s(x2) + s(x3),
  data = d, family = "normal")
plot(b)

## Compute the log-likelihood using the family object.
f <- family(b)
sum(f$d(y = d$num, par = f$map2par(fitted(b)), log = TRUE))

## For using JAGS() more details are needed.
norm4JAGS.bamlss <- function(...) {
  f <- normal.bamlss()
  f$bugs <- list(
    "dist" = "dnorm",
    "eta" = BUGSeta,
    "model" = BUGSmodel,
    "reparam" = c(sigma = "1 / sqrt(sigma)")
  )
  return(f)
}

## Now with bfit() and JAGS().
b <- bamlss(num ~ s(x1) + s(x2) + s(x3), data = d,
  optimizer = bfit, sampler = JAGS, family = "norm4JAGS")
plot(b)
}
}

\keyword{regression}

