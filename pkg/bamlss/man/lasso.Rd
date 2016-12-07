\name{la}
\alias{la}
\alias{lasso}
\alias{lasso.plot}
\alias{lasso.stop}

\title{Lasso Smooth Constructor}

\description{
  Smooth constructors and optimizer for Lasso penalization with \code{\link{bamlss}}. The
  penalization is based on a Taylor series approximation of the Lasso penalty.
}

\usage{
## Smooth constructor function.
la(formula, type = c("single", "multiple"), ...)

## Single Lasso smoothing parameter optimizer.
lasso(x, y, start = NULL, lower = 0.001, upper = 1000,
  nlambda = 100, lambda = NULL, verbose = TRUE, digits = 4,
  flush = TRUE, ...)

## Plotting function for lasso() optimizer.
lasso.plot(x, which = c("criterion", "parameters"), spar = TRUE, ...)

## Extract optimum stopping iteration for lasso() optimizer.
## Based on the minimum of the information criterion.
lasso.stop(x)
}

\arguments{
  \item{formula}{A formula like \code{~ x1 + x2 + ... + xk} of variables which should be
    penalized with Lasso.}
  \item{type}{Should one single penalty parameter be used or multiple parameters, one for each
    covariate in \code{formula}.}
  \item{x}{For function \code{lasso()} the \code{x} list, as returned from function
    \code{\link{bamlss.frame}}, holding all model matrices and other information that is used for
    fitting the model. For the plotting function and \code{lasso.stop()} the corresponding
    \code{\link{bamlss}} object fitted with the \code{lasso()} optimizer.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{start}{A vector of starting values. Note, Lasso smoothing parameters will be dropped.}
  \item{lower}{Numeric. The minimum lambda value.}
  \item{upper}{Numeric. The maximum lambda value.}
  \item{nlambda}{Integer. The number of smoothing parameters for which coefficients should be
    estimated, i.e., the vector of smoothing parameters is build up as a sequence from
    \code{lower} to \code{upper} with length \code{nlambda}.}
  \item{lambda}{Numeric. A sequence/vector of lambda parameters that should be used.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{digits}{Set the digits for printing when \code{verbose = TRUE}.}
  \item{flush}{use \code{\link{flush.console}} for displaying the current output in the console.}
  \item{which}{Which of the two provided plots should be created?}
  \item{spar}{Should graphical parameters be set by the plotting function?}
  \item{\dots}{Arguments passed to the subsequent smooth constructor function.
    \code{lambda} controls the starting value of the penalty parameter, \code{const} the constant
    that is added within the penalty approximation.
    For the optimizer \code{lasso()} arguments passed to function \code{\link{bfit}}.}
}

\value{
  For function \code{la()}, similar to function \code{\link[mgcv]{s}} a simple smooth
  specification object.

  For function \code{lasso()} a list containing the following objects:
  \item{fitted.values}{A named list of the fitted values based on the last lasso iteration
    of the modeled parameters of the selected distribution.}
  \item{parameters}{A matrix, each row corresponds to the parameter values of one boosting iteration.}
  \item{lasso.stats}{A matrix containing information about the log-likelihood, log-posterior
    and the information criterion for each lambda.}
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
\dontrun{## Load data example.
data(rent99, package = "gamlss.data")
rent99$district <- as.factor(rent99$district)

## Model formula.
fl <- ~ area + yearc + location + bath + kitchen + cheating + district
f <- list(rentsqm ~ la(fl), sigma ~ la(fl))

## Estimate Lasso model with
## single lambda parameter.
b1 <- bamlss(f, data = rent99, sampler = FALSE, optimizer = lasso,
  lower = 0.01, upper = 1, criterion = "BIC", nlambda = 100)

## Plot information criterion and
## coefficient paths.
lasso.plot(b1)

## Extract optimum lambda/iteration.
lasso.stop(b1)

## Extract coefficients for optimum
## Lasso parameter.
parameters(b1, mstop = lasso.stop(b1))

## Predict with optimum lasso.
p1 <- predict(b1, mstop = lasso.stop(b1))

## Now with different lambda for mu and sigma using
## optimizer function bfit().
b2 <- bamlss(f, data = rent99, criterion = "BIC", sampler = FALSE)
summary(b2)

## With lambda parameters for each covariate.
f <- list(
  rentsqm ~ la(area) + la(yearc) + la(location) +
    la(bath) + la(kitchen) + la(cheating) + la(district),
  sigma ~ la(area) + la(yearc) + la(location) +
    la(bath) + la(kitchen) + la(cheating) + la(district)
)
b3 <- bamlss(f, data = rent99, criterion = "BIC", sampler = FALSE)
summary(b3)

## With full MCMC.
b4 <- bamlss(f, data = rent99, criterion = "BIC")
summary(b4)
}
}

\keyword{regression}

