\name{bfit}
\alias{bfit}

\title{
  Fit a BAMLSS with Backfitting
}

\description{
  This optimizer function is a generic tool for fitting BAMLSS using a backfitting algorithm.
  The backfitting procedure is based on iteratively reweighted least squares (IWLS), however,
  the updating methods for model terms can be more general, see the details. In addition, the
  default IWLS updating scheme implements optimum smoothing variance selection based on information
  criteria using a stepwise approach.
}

\usage{
bfit(x, y, family, start = NULL, weights = NULL, offset = NULL,
  update = "iwls", criterion = c("AICc", "BIC", "AIC"),
  eps = .Machine$double.eps^0.25, maxit = 400,
  outer = FALSE, inner = FALSE, mgcv = FALSE,
  verbose = TRUE, digits = 4, ...)
}

\arguments{
  \item{x}{The \code{x} list, as returned from function \code{\link{bamlss.frame}}, holding all model
    matrices and other information that is used for fitting the model.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{family}{A \pkg{bamlss} family object, see \code{\link{family.bamlss}}.}
  \item{start}{A named numeric vector containing possible starting values, the names are based on
    function \code{\link{parameters}}.}
  \item{weights}{Prior weights on the data, as returned from function \code{\link{bamlss.frame}}.}
  \item{offset}{Can be used to supply model offsets for use in fitting,
    returned from function \code{\link{bamlss.frame}}.}
  \item{update}{Sets the updating function for model terms, per default \code{"iwls"}, a
    character pointing to the internal set of updating functions. Other options are
    \code{"optim"} and \code{"lm"}, however this is more experimental
    and should not be set by the user. Another option is to pass an updating function which
    should be used for each model term, the structure of updating functions is described in the
    details below.}
  \item{criterion}{Set the information criterion that should be used, e.g., for smoothing
    variance selection.}
  \item{eps}{The relative convergence tolerance of the backfitting algorithm.}
  \item{maxit}{The maximum number of iterations for the backfitting algorithm}
  \item{outer}{Should the current working observations be computed in one outer iteration,
    otherwise the working observations are computed anew for each model term updating step.}
  \item{inner}{Should the model terms for one parameter of the modeled distribution be
    fully updated until convergence in an inner iteration, i.e., the algorithm waits until
    coefficients for the current distribution parameter do not change anymore before updating
    the next parameter.}
  \item{mgcv}{Should the \pkg{mgcv} \code{\link{gam}} function be used for computing updates in
    an \code{inner} iteration with working observations provided in an \code{outer} iteration.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{digits}{Set the digits for printing when \code{verbose = TRUE}.}
  \item{\dots}{Arguments passed to function \code{\link{bamlss.engine.setup}}.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
  A list containing the following objects:
  \item{fitted.values}{A named list of the fitted values of the modeled parameters
    of the selected distribution.}
  \item{parameters}{The estimated set regression coefficients and smoothing variances.}
  \item{edf}{The equivalent degrees of freedom used to fit the model.}
  \item{logLik}{The value of the log-likelihood.}
  \item{logPost}{The value of the log-posterior.}
  \item{IC}{The value of the information criterion.}
  \item{converged}{Logical, indicating convergence of the backfitting algorithm.}
}

\references{
  Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location,
  scale and shape, (with discussion), \emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

  Belitz C, Lang S (2008). Simultaneous Selection of Variables and Smoothing Parameters in
  Structured Additive Regression Models. \emph{Computational Statistics \& Data Analysis},
  \bold{53}, pp 61-81.
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}},
  \code{\link{bamlss.engine.setup}}, \code{\link{set.starting.values}}
}

\examples{
## Simulated data example illustrating
## how to call the optimizer function.
## This is done internally within
## the setup of function bamlss().
d <- GAMart()
f <- num ~ s(x1) + s(x2) + s(x3) + te(lon, lat)
bf <- bamlss.frame(f, data = d, family = "gaussian")
opt <- with(bf, bfit(x, y, family))
print(str(opt))
}

\keyword{regression}

