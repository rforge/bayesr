\name{jm.mode}
\alias{jm.mode}

\title{
  Joint model posterior mode estimation
}

\description{
  This function computes posterior mode estimates of the parameters of a flexible joint model
  of longitudinal and survival data with structured additive predictors using a Newton-Raphson 
  algorithm. Integrals are solved numerically using trapezoidal rule. Smoothing variances are computed using numerical
  optimization.
}

\usage{
jm.mode(x, y, weights, offset,
  criterion = c("AICc", "BIC", "AIC"), maxit = c(100, 1),
  nu = c("lambda" = 0.1, "gamma" = 0.1, "mu" = 0.1, "sigma" = 0.1, "alpha" = 0.1, "dalpha" = 0.1),
  update.nu = TRUE,  
  eps = 0.0001, nu.eps = -Inf, edf.eps = NULL, ic.eps = 1e-08, alpha.eps = 0.001, 
  nback = 40, verbose = TRUE, digits = 4, ...)
}

\arguments{
 \item{x}{The \code{x} list, as returned from function
    \code{\link{bamlss.frame}} and transformed by function \code{\link{jm.transform}},
    holding all model matrices and other information that is used for
    fitting the model.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{weights}{Prior weights on the data, as returned from function \code{\link{bamlss.frame}}.}
  \item{offset}{Can be used to supply model offsets for use in fitting,
    returned from function \code{\link{bamlss.frame}}.}
  \item{criterion}{Information criterion to be used, e.g., for smoothing
    variance selection. Options are the corrected AIC \code{"AICc"} (see Details), the \code{"BIC"} and
    \code{"AIC"}. Defaults to \code{"AICc"}?}
  \item{maxit}{Vector containing the maximum number of iterations for the backfitting algorithm with \code{maxit[1]} defining the iterations for the full model and \code{maxit[2]} the iterations within each predictor. \code{maxit[2]} defaults to 1 if only one value is specified.}
  \item{nu}{Vector of step lengths for parameter updates of one Newton-Raphson update for 
  each predictor of the joint model.}
  \item{update.nu}{Should the updating step length be optimized in each iteration
    of the backfitting algorithm? Overridess parameter \code{nu} if set to \code{TRUE}.}
  \item{eps}{The relative convergence tolerance of the backfitting algorithm.}
  \item{nu.eps}{The relative convergence tolerance for the optimization of steplength \code{nu}.}
  \item{ic.eps}{Raus?}
  \item{alpha.eps}{The relative convergence tolerance of the backfitting algorithm for predictor alpha.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{digits}{Set the digits for printing when \code{verbose = TRUE}.}
  \item{\dots}{Currently not used.}
}

\value{
  A list containing the following objects:
  \item{fitted.values}{A named list of the fitted values of the modelled parameters
    of the selected distribution.}
  \item{parameters}{The estimated set of regression coefficients and smoothing variances.}
  \item{edf}{The equivalent degrees of freedom used to fit the model.}
  \item{logLik}{The value of the log-likelihood.}
  \item{logPost}{The value of the log-posterior.}
  \item{hessian}{The Hessian matrix evaluated at the posterior mode.}
  \item{converged}{Logical, indicating convergence of the backfitting algorithm.}
  \item{time}{The runtime of the algorithm.}
  \item{intercept_paths}{Raus?}
  \item{eta_timegrid}{Raus?}
}

\references{
  Umlauf N, Klein N, Zeileis A (2016). Bayesian Additive Models for Location
    Scale and Shape (and Beyond). \emph{(to appear)}
}

\details{
(1) Numerical optimization for AIC and nu.
(2) Long und Surv erstmal allein optimiert. 
(3) AICc Hurvich
}

\seealso{
  \code{\link{jm.mcmc}}, \code{\link{jm.bamlss}}, \code{\link{jm.transform}},
  \code{\link{jm.sim}}, \code{\link{bamlss}}
}

\examples{
## Please see the examples of function jm.mcmc()!
}

\keyword{regression}
\keyword{survival}

