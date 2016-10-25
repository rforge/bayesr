\name{boost}
\alias{boost}
\alias{boost.summary}
\alias{boost.summary}
\alias{print.boost.summary}
\alias{plot.boost.summary}

\title{Gradient Boosting BAMLSS}

\description{
  Optimizer function for gradient boosting with \code{\link{bamlss}}. In each boosting iteration
  the function selects the model term with the largest contribution to the log-likelihood.
}

\usage{
## Gradient boosting optimizer:
boost(x, y, family,
  nu = 0.1, df = 4, maxit = 100, mstop = NULL,
  verbose = TRUE, digits = 4, flush = TRUE,
  nback = 10, Machine$double.eps^0.25,
  plot = TRUE, ...)

## Boosting summary extractor:
boost.summary(object, ...)

## Boosting summary printing and plotting:
\method{print}{boost.summary}(x, summary = TRUE, plot = TRUE, ...)
\method{plot}{boost.summary}(x, ...)
}

\arguments{
  \item{x}{For function \code{boost()} the \code{x} list, as returned from function
    \code{\link{bamlss.frame}}, holding all model matrices and other information that is used for
    fitting the model. For the plotting function the corresponding \code{\link{bamlss}} object
    fitted with the \code{boost()} optimizer.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{family}{A \pkg{bamlss} family object, see \code{\link{family.bamlss}}.}
  \item{nu}{Numeric, between [0, 1], controls the step size, i.e., the amount
    that should be added to model term parameters.}
  \item{df}{Integer, defines the initial degrees of freedom that should be assigned
    to each smooth model term. May also be a named vector, the names must match the model term
    labels, e.g., as provided in \code{\link{summary.bamlss}}.}
  \item{maxit}{Integer, the maximum number of boosting iterations.}
  \item{mstop}{For convenience, overwrites \code{maxit}.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{digits}{Set the digits for printing when \code{verbose = TRUE}.}
  \item{flush}{use \code{\link{flush.console}} for displaying the current output in the console.}
  \item{nback}{If \code{nback} is not \code{NULL}, then the algorithm stops if the difference
    of relative changes of predictors of the last \code{nback} iterations is smaller or equal to
    \code{eps}. If argument \code{mstop} is \code{NULL} the default is to set relatively large to
    \code{mstop = 10000}.}
  \item{eps}{The relative convergence tolerance of the boosting algorithm using \code{nback}.}
  \item{plot}{Should the boosting summary be printed and plotted?}
  \item{object}{A \code{\link{bamlss}} object that was fitted using \code{boost()}.}
  \item{summary}{Should the summary be printed?}
  \item{\dots}{For function \code{boost()}, arguments passed to \code{\link{bamlss.engine.setup}}.
    for function \code{boost.summary()} arguments passed to function \code{print.boost.summary()}.}
}

\value{
  For function \code{boost.summary()} a list containing information on selection frequencies etc.
  For function \code{boost()} a list containing the following objects:
  \item{fitted.values}{A named list of the fitted values based on the last boosting iteration
    of the modeled parameters of the selected distribution.}
  \item{parameters}{A matrix, each row corresponds to the parameter values of one boosting iteration.}
  \item{boost.summary}{The boosting summary which can be printed and plotted.}
}

\section{WARNINGS}{
  The function is currently experimental!
  The function does not take care of variable scaling! This must be done by the user, e.g.,
  one option is to use argument \code{scale.d} in function \code{\link{bamlss.frame}}.
}

\seealso{
\code{\link{bamlss.frame}}, \code{\link{bamlss}}
}

\examples{
\dontrun{## Simulate data.
set.seed(123)
d <- GAMart()

## Estimate model.
f <- num ~ x1 + x2 + x3 + lon + lat +
  s(x1) + s(x2) + s(x3) + s(lon) + s(lat) + te(lon,lat)
b <- bamlss(f, data = d, optimizer = boost,
  sampler = FALSE, scale.d = TRUE)

## Plot estimated effects.
plot(b)

## Print and plot the boosting summary.
boost.summary(b)

## Extract estimated parameters for certain
## boosting iterations.
parameters(b, mstop = 1)
parameters(b, mstop = 100)

## Also works with predict().
head(do.call("cbind", predict(b, mstop = 1)))
head(do.call("cbind", predict(b, mstop = 100)))
}
}

\keyword{regression}

