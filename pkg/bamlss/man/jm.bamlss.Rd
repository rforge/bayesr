\name{jm.bamlss}
\alias{jm.mode}
\alias{jm.mcmc}


\title{
  Fit flexible additive joint models in \pkg{bamlss}
}

\description{
  Family object to fit a flexible additive joint model for longitudinal and survival
  data under a Bayesian approach as presented in Koehler, Umlauf and Greven (2016). 
  All parts of the joint model can be specified as structured additive predictors. 
  See the details and examples.
}

\usage{
bamlss(formula, family = "jm", data = NULL, timevar=NULL,
  idvar=NULL, start = NULL, knots = NULL, weights = NULL,
  subset = NULL, offset = NULL, na.action = na.omit,
  contrasts = NULL, reference = NULL, transform = jm.transform,
  optimizer = jm.mode, sampler = jm.mcmc, samplestats = NULL,
  results = NULL, cores = NULL, sleep = NULL,
  combine = TRUE, model = TRUE, x = TRUE, ...)
  
jm.bamlss(...)
}

\arguments{
  \item{formula}{list of formulas where each list entry specifies the details of one predictor in
    the joint model. For incorporating smooth terms, all model term constructors
    implemented in \code{\link{mgcv}} such as \code{\link[mgcv]{s}}, 
    \code{\link[mgcv]{te}} and \code{\link[mgcv]{ti}} can be used, amongst others. 
    See Examples.}
  \item{family}{A \code{\link{bamlss.family}} object, specifying the details of the modeled
     distribution such as the parameter names, the density function, link functions, 
     etc. \code{"jm"} specifies joint model.}
  \item{data}{A \code{\link{data.frame}} or \code{\link{list}} containing the model
    response variable(s) and covariates specified in the \code{formula} in long format.
    By default the variables are taken from \code{environment(formula)}:
    typically the environment from which \code{bamlss} is called.}
  \item{timevar}{name of the time variable in the longitudinal model part}
  \item{idvar}{name of the subject's identifier connecting multiple measurements of the same
    subject}
  \item{start}{A named numeric vector containing starting values.  If 
    \code{NULL} posterior mode estimates are used as starting values.}
  \item{knots}{An optional list containing user specified knots, see the documentation of
    function \code{\link[mgcv]{gam}}.}
  \item{weights}{Prior weights on the data.}
  \item{subset}{An optional vector specifying a subset of observations to be
    used in the fitting process.}
  \item{offset}{Can be used to supply model offsets for use in fitting.}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}'s. The default is set by the \code{na.action}
    setting of \code{\link{options}}, and is \code{\link{na.omit}} if set to \code{NULL}.}
  \item{contrasts}{An optional list. See the \code{contrasts.arg} of 
    \code{\link[stats]{model.matrix.default}}.}
  \item{reference}{A \code{character} specifying a reference category, e.g., when
    fitting a multinomial model.}
  \item{transform}{A transformer function that is applied on the \code{\link{bamlss.frame}}.
    See, e.g., function \code{\link{randomize}} and \code{\link{bamlss.engine.setup}}.}
  \item{optimizer}{An optimizer function that returns, e.g., posterior mode estimates
    of the parameters as a named numeric vector. The default optimizer function is
    \code{\link{bfit}}. If set to \code{FALSE}, no optimizer function will be used.}
  \item{sampler}{A sampler function that returns a matrix of samples, the columns represent the
    parameters, the rows the iterations. The returned matrix must be coerced to an object of
    class \code{"mcmc"}, see \code{\link[coda]{as.mcmc}}. The default sampler function is
    \code{\link{GMCMC}}. If set to \code{FALSE}, no sampler function will be used.}
  \item{samplestats}{A function computing statistics from samples, per default function
    \code{\link{samplestats}} is used. If set to \code{FALSE}, no \code{samplestats} function
    will be used. Note that this option is crucial for very large datasets, as computing
    statistics from samples this way may be very time consuming!}
  \item{results}{A function computing results from the parameters and/or samples, e.g., for
    creating effect plots, see function \code{link{results.bamlss.default}}. If set \code{FALSE}
    no \code{results} function will be used.}
  \item{cores}{An integer specifying the number of cores that should be used for the sampler
    function. This is based on function \code{\link[parallel]{mclapply}} of the \pkg{parallel}
    package.}
  \item{sleep}{Time the system should sleep before the next core is started.}
  \item{combine}{If samples are computed on multiple cores, should the samples be combined into
    one bin \code{\link[coda]{mcmc}} matrix?}
  \item{model}{If set to \code{FALSE} the model frame used for modeling is not part of the
    return value.}
  \item{x}{If set to \code{FALSE} the model matrices are not part of the return value.}
  \item{\dots}{Arguments passed to the \code{transformer}, \code{optimizer}, \code{sampler},
    \code{results} and \code{samplestats} function.}
}



\details{
Survival model

Longitudinal model

Each predictor is sum of functions

Integration using trapezoidal rule

Posterior mode estimation

Posterior mean estimation


}


\references{
  Koehler N, Umlauf N, Greven S (2016). Flexible Bayesian additive joint models with an application     to type 1 diabetes research. \emph{(submitted)}
  
  Umlauf N, Klein N, Zeileis A (2016). Bayesian Additive Models for Location
    Scale and Shape (and Beyond). \emph{(to appear)}
}


\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}}
}

\examples{
set.seed(123)

## Simulate survival data
# with random intercepts and linear effect of time,
# constant association alpha and no effect of the derivative
d <- simJM(nsub = 200, long_setting = "linear", alpha_setting = "constant",
            dalpha_setting = "zero", full=FALSE)

## Formula of the according joint model
f <- list(Surv2(survtime, event, obs = y) ~ s(survtime, bs="ps"),
          gamma ~ s(x1, bs="ps"),
          mu ~ obstime + s("id", bs="re"),
          sigma ~ 1,
          alpha ~ 1,
          dalpha ~ -1)

## Joint model estimation
# The family object jm.bamlss sets the default optimizer and sampler function.
# First, posterior mode estimates are computed using function jm.mode(), 
# afterwards the sampler jm.mcmc() is started.
b <- bamlss(f, data=d, family="jm", timevar="obstime", idvar="id")

## Plot estimated effects.
plot(b)



## Simulate survival data
# with functional random intercepts and a nonlinear effect of time,
# time-varying association alpha and no effect of the derivative.
# Note: This specification is the default in the data simulation.
d <- simJM(nsub = 200, full=FALSE)

## Formula of the according joint model
# specify the number of knots for the smooth nonlinear effect of time
klong <- 8
f <- list(Surv2(survtime, event, obs = y) ~ s(survtime, bs="ps"),
          gamma ~ s(x1, bs="ps"),
          mu ~ ti(id, bs="re") + ti(obstime, bs="ps", k=klong) + 
               ti(id, obstime, bs=c("re", "ps"), k=c(nlevels(d$id), klong)) + s(x2, bs="ps"),
          sigma ~ 1,
          alpha ~ s(survtime, bs="ps"),
          dalpha ~ -1)

## Estimating posterior mode only using jm.mode()
b_mode <- bamlss(f, data=d, family="jm", timevar="obstime", idvar="id", 
             sampler = MVNORM)

## Estimating posterior means using jm.mcmc() and using mode starting values
b_mean <- bamlss(f, data=d, family="jm", timevar="obstime", idvar="id",
                 optimizer=FALSE, jm.start = parameters(b_mode))

}

\keyword{regression}

