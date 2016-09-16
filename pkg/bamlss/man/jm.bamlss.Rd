\name{jm_bamlss}
\alias{jm_bamlss}
\alias{jm.mode}
\alias{jm.mcmc}
\alias{jm.transform}

\title{
  Fit Flexible Additive Joint Models
}

\description{
  Family object to fit a flexible additive joint model for longitudinal and survival
  data under a Bayesian approach as presented in Koehler, Umlauf and Greven (2016). 
  All parts of the joint model can be specified as structured additive predictors. 
  See the details and examples.
}

\usage{
## JM family object.
jm_bamlss(...)

## "bamlss.frame" transformer function
## to set up joint models.
jm.transform(x, y, data, terms, knots, formula,
  family, subdivisions = 25,
  timedependent = c("lambda", "mu", "alpha", "dalpha"),
  timevar = NULL, idvar = NULL, alpha = .Machine$double.eps,
  mu = NULL, sigma = NULL, sparse = TRUE, ...)

## Posterior mode optimizing engine.     
jm.mode(x, y, start = NULL, weights = NULL,
  offset = NULL, criterion = c("AICc", "BIC", "AIC"),
  maxit = c(100, 1), nu = c("lambda" = 0.1, "gamma" = 0.1,
    "mu" = 0.1, "sigma" = 0.1, "alpha" = 0.1, "dalpha" = 0.1),
  update.nu = TRUE, eps = 0.0001, alpha.eps = 0.001,
  ic.eps = 1e-08, nback = 40, verbose = TRUE,
  digits = 4, ...)

## Sampler function.
jm.mcmc(x, y, family, start = NULL, weights = NULL,
  offset = NULL, n.iter = 1200, burnin = 200, thin = 1,
  verbose = TRUE, digits = 4, step = 20, ...)
}

\arguments{
  \item{x}{The \code{x} list, as returned from function
    \code{\link{bamlss.frame}} (and transformed by function \code{jm.transform()}),
    holding all model matrices and other information that is used for
    fitting the model.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{data}{A \code{\link{data.frame}} or \code{\link{list}} containing the model
    response variable(s) and covariates specified in the \code{formula} in long format.
    By default the variables are taken from \code{environment(formula)}:
    typically the environment from which \code{bamlss} is called.}
  \item{terms}{The corresponding \code{\link{terms.bamlss}} object needed for processing.}
  \item{knots}{An optional list containing user specified knots, see the documentation of
    function \code{\link[mgcv]{gam}}.}
  \item{formula}{The corresponding \code{\link{bamlss.formula}}.}
  \item{family}{The \code{\link{bamlss.family}} object.}
  \item{subdivisions}{How many time points should be created for each individual.}
  \item{timedependent}{A character vector specifying the names of parameters in \code{x}
    that are time-dependent. Time grid design matrices are only computed for these parameters.}
  \item{timevar}{A character specifying the name of the survival time variable in the
    data set.}
  \item{idvar}{Depending on the type of data set, this is the name of the variable specifying
     identifier of individuals.}
  \item{alpha}{Numeric, a starting value for the intercept of the association parameter alpha.}
  \item{mu}{Numeric, a starting value for the intercept of the mu parameter.}
  \item{sigma}{Numeric, a starting value for the intercept of the sigma parameter.}
  \item{sparse}{Logical, should sparse matrix structures be used for updating and sampling
    of mu parameter model terms?}
  \item{start}{A named numeric vector containing possible starting values, the names are based on
    function \code{\link{parameters}}.}
  \item{weights}{Currently not supported.}
  \item{offset}{Currently not supported.}
  \item{criterion}{Information criterion to be used, e.g., for smoothing
    variance selection. Options are the corrected AIC \code{"AICc"} (see Details), the \code{"BIC"} and
    \code{"AIC"}. Defaults to \code{"AICc"}?}
  \item{maxit}{Vector containing the maximum number of iterations for the backfitting
    algorithm with \code{maxit[1]} defining the iterations for the full model and \code{maxit[2]}
    the iterations within each predictor. \code{maxit[2]} defaults to 1 if only one value is
    specified.}
  \item{nu}{Vector of step lengths for parameter updates of one Newton-Raphson update for 
    each predictor of the joint model.}
  \item{update.nu}{Should the updating step length be optimized in each iteration
    of the backfitting algorithm? Uses \code{nu} as starting value if set to \code{TRUE}.}
  \item{eps}{The relative convergence tolerance of the backfitting algorithm.}
  \item{alpha.eps}{The relative convergence tolerance of the backfitting algorithm
    for predictor alpha.}
  \item{ic.eps}{The relative convergence tolerance of the information criterion used, e.g.,
    for smoothing variance selection.}
  \item{nback}{For computing \code{ic.eps}, how many iterations back should be included
    when computing relative convergence tolerance of the information criterion.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{digits}{Set the digits for printing when \code{verbose = TRUE}.}
  \item{n.iter}{the number of MCMC iterations.}
  \item{burnin}{the burn-in phase of the sampler, i.e., the number of starting samples that
    should be removed.}
  \item{thin}{the thinning parameter for MCMC simulation. E.g., \code{thin = 10} means,
    that only every 10th sampled parameter will be stored.}
  \item{step}{How many times should algorithm runtime information be printed, divides \code{n.iter}.}
  \item{\dots}{Currently not used.}
}

\details{
We refer to the paper of Koehler, Umlauf, and Greven (2016) for details on the flexible
additive joint model. In short, we model the hazard of subject \eqn{i} an event at time
\eqn{t} as
\deqn{h_{i}(t) = \exp [\eta_{\lambda i}(t) + \eta_{\gamma i}+\eta_{\alpha i}(t) \cdot \eta_{\mu i}(t)]}
with predictor \eqn{\eta_{\lambda}} for all survival covariates that are time-varying or have a
time-varying coefficient (including the log baseline hazard), predictor \eqn{\eta_{\gamma}} for
baseline survival covariates, predictor \eqn{\eta_{\alpha}} representing the potentially
time-varying association between the longitudinal marker \eqn{\eta_{\mu}} and the hazard. The
longitudinal response \eqn{y_{ij}} at time points \eqn{t_{ij}} is modeled as
\deqn{y_{ij}=\eta_{\mu i}(t_{ij})+e_{ij}} with independent normal errors \eqn{N(0, \exp[\eta_{\sigma i}(t_{ij})]^2)}. 

Each predictor \eqn{\eta_{ki}} is a structured additive predictor, i.e. a sum of functions of
covariates \eqn{\eta_{ki} =  \sum_{m=1}^{M_k} f_{km}(\bm{x}_{ki})}. Each of these functions can be
modeled parametrically or using basis function evaluations from the smooth constructors in
\pkg{mgcv} such as \code{\link[mgcv]{s}}, \code{\link[mgcv]{te}} and \code{\link[mgcv]{ti}} and can
include smooth time-varying, random or spatial effects. For the Bayesian estimation of these effects
we specify corresponding priors: For linear or parametric terms we use vague normal priors, smooth
and random effect terms are regularized by placing generic multivariate normal priors on the
coefficients and for anisotropic smooths, when multiple smoothing variance parameters are involved,
more complex prior are in place (cf. Koehler, Umlauf & Greven, 2016). We use inverse Gamma
hyper-priors, i.e. IG(0.001, 0.001) to obtain an inverse Gamma full conditional for the variance
parameters. We estimate the posterior mode by maximizing the log-posterior of the model using a
Newton-Raphson procedure, the posterior mean is obtained via derivative-based Metropolis-Hastings
sampling. We recommend to use posterior mode estimates for a quick model assessment. In order to
draw correct inferences from the model, posterior mean estimates should be computed. 
We approximate integration in the survival part of the likelihood using trapezoidal rule. For
posterior mode estimation.
}

\references{
  Koehler N, Umlauf N, Greven S (2016). Flexible Bayesian Additive Joint Models with an Application
    to Type 1 Diabetes Research. \emph{(submitted)}
  
  Umlauf N, Klein N, Zeileis A (2016). Bayesian Additive Models for Location
    Scale and Shape (and Beyond). \emph{(to appear)}
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}}.
}

\examples{
\dontrun{set.seed(123)

## Simulate survival data
## with random intercepts/slopes and a linear effect of time,
## constant association alpha and no effect of the derivative
d <- simJM(nsub = 200, long_setting = "linear", 
  alpha_setting = "constant",
  dalpha_setting = "zero", full = FALSE)

## Formula of the according joint model
f <- list(
  Surv2(survtime, event, obs = y) ~ s(survtime, bs = "ps"),
  gamma ~ s(x1, bs = "ps"),
  mu ~ obstime + s(id, bs = "re") + 
    s(id, obstime, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

## Joint model estimation
## jm_bamlss() sets the default optimizer and sampler function.
## First, posterior mode estimates are computed using function
## jm.mode(), afterwards the sampler jm.mcmc() is started.
b <- bamlss(f, data = d, family = "jm",
  timevar = "obstime", idvar = "id")

## Plot estimated effects.
plot(b)

## Predict event probabilities for two
## individuals.
nd <- subset(d, id %in% c(192, 127))
p <- predict(b, newdata = nd, type = "probabilities")
print(p)

## Plot of survival probabilities and
## corresponding longitudinal effects
## for individual i.
plot_jm_surv <- function(i, maxtime = 80, grid = 10) {
  on.exit(par(par(no.readonly = TRUE)))

  nd <- data.frame(
    "id" = factor(i, levels = levels(d$id)),
    "x1" = d$x1[d$id == i][1]
  )
  time <- seq(0, maxtime, length = grid)

  if(all(time == time[1]))
    stop(paste("Not enough time points available for individual ", i, "!", sep = ""))

  p_surv <- NULL
  for(t in time) {
    nd$survtime <- t
    nd$obstime <- t
    p_surv <- rbind(p_surv, t(predict(b, newdata = nd, type = "probabilities", FUN = c95)))
  }

  nd <- data.frame(
    "obstime" = time,
    "id" = factor(i, levels = levels(d$id))
  )
  p_long <- predict(b, newdata = nd, model = "mu", FUN = c95)

  par(mfrow = c(2, 1), mar = rep(0, 4), oma = c(4.1, 4.1, 1.1, 4.1))
  plot2d(p_surv ~ time, fill.select = c(0, 1, 0, 1),
    scheme = 2, axes = FALSE, ylim = c(0, 1))
  axis(2)
  box()
  mtext("Prob(T > t)", side = 2, line = 2.5)
  plot2d(p_long ~ time, fill.select = c(0, 1, 0, 1), scheme = 2,
    xlab = "Time", ylab = "Effect of time", axes = FALSE)
  with(subset(d, id == i), points(obstime, y))
  axis(1)
  axis(4)
  box()
  mtext("Effect of time", side = 4, line = 2.5)
  mtext("Time", side = 1, line = 2.5)
}

plot_jm_surv(i = 3)
plot_jm_surv(i = 30)
plot_jm_surv(i = 192)

## Simulate survival data
## with functional random intercepts and a nonlinear effect 
## of time, time-varying association alpha and no effect 
## of the derivative.
## Note: This specification is the simJM default.
d <- simJM(nsub = 200, full = FALSE)

## Formula of the according joint model
## umber of knots for the smooth nonlinear effect of time
klong <- 8
f <- list(
  Surv2(survtime, event, obs = y) ~ s(survtime, bs = "ps"),
  gamma ~ s(x1, bs = "ps"),
  mu ~ ti(id, bs = "re") + 
    ti(obstime, bs = "ps", k = klong) + 
    ti(id, obstime, bs = c("re", "ps"), 
      k = c(nlevels(d$id), klong)) + 
    s(x2, bs = "ps"),
  sigma ~ 1,
  alpha ~ s(survtime, bs = "ps"),
  dalpha ~ -1
)

## Estimating posterior mode only using jm.mode()
b_mode <- bamlss(f, data = d, family = "jm",
  timevar = "obstime", idvar = "id", sampler = MVNORM)

## Estimating posterior means using jm.mcmc() 
##  with starting values generated from posterior mode 
b_mean <- bamlss(f, data = d, family = "jm",
  timevar = "obstime", idvar = "id", optimizer = FALSE,
  jm.start = parameters(b_mode))

## Plot effects.
plot(b_mean, model = "alpha")
plot(b_mean, model = "alpha", which = samples)
}
}

\keyword{regression}
\keyword{models}

