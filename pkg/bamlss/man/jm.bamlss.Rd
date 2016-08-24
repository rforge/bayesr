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

jm.bamlss(formula, data = NULL, timevar=NULL, idvar = NULL, 
          start = NULL, knots = NULL, weights = NULL,
          offset = NULL, cores = NULL, combine = TRUE, ...)
          
jm.mode(x, y, weights, offset, criterion = c("AICc", "BIC", "AIC"), 
        maxit = c(100, 1), nu = c("lambda" = 0.1, "gamma" = 0.1, 
        "mu" = 0.1, "sigma" = 0.1, "alpha" = 0.1, "dalpha" = 0.1),
        update.nu = TRUE, eps = 0.0001, alpha.eps = 0.001,
        ic.eps = 1e-08, nback = 40, verbose = TRUE,
        digits = 4, ...)

jm.mcmc(x, y, family, start, jm.start = NULL,
        weights, offset, n.iter = 1200, burnin = 200,
        thin = 1, verbose = TRUE, digits = 4,
        step = 20, ...)
}


\arguments{
  \item{formula}{list of formulas where each list entry specifies the details of one predictor in
    the joint model. For incorporating smooth terms, all model term constructors
    implemented in \pkg{mgcv} such as \code{\link[mgcv]{s}}, 
    \code{\link[mgcv]{te}} and \code{\link[mgcv]{ti}} can be used, amongst others. 
    See Examples.}
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
  \item{offset}{Can be used to supply model offsets for use in fitting.}
  \item{cores}{An integer specifying the number of cores that should be used for the sampler
    function. This is based on function \code{\link[parallel]{mclapply}} of the \pkg{parallel}
    package.}
  \item{combine}{If samples are computed on multiple cores, should the samples be combined into
    one bin \code{\link[coda]{mcmc}} matrix?}
  \item{\dots}{Arguments passed to the functions as \code{jm.transform} and \code{jm.mode}, \code{jm.mcmc}. For further options in the bamlss call refer to \code{\link{bamlss}}.}
  \item{x}{The \code{x} list, as returned from function
    \code{\link{bamlss.frame}} and transformed by function \code{jm.transform},
    holding all model matrices and other information that is used for
    fitting the model.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{criterion}{Information criterion to be used, e.g., for smoothing
    variance selection. Options are the corrected AIC \code{"AICc"} (see Details), the \code{"BIC"} and
    \code{"AIC"}. Defaults to \code{"AICc"}?}
  \item{maxit}{Vector containing the maximum number of iterations for the backfitting algorithm with \code{maxit[1]} defining the iterations for the full model and \code{maxit[2]} the iterations within each predictor. \code{maxit[2]} defaults to 1 if only one value is specified.}
  \item{nu}{Vector of step lengths for parameter updates of one Newton-Raphson update for 
  each predictor of the joint model.}
  \item{update.nu}{Should the updating step length be optimized in each iteration
    of the backfitting algorithm? Uses \code{nu} as starting value if set to \code{TRUE}.}
  \item{eps}{The relative convergence tolerance of the backfitting algorithm.}
  \item{alpha.eps}{The relative convergence tolerance of the backfitting algorithm for predictor alpha.}
  \item{ic.eps}{The relative convergence tolerance of the information criterion used, e.g. for smoothing variance selection. }
  \item{nback}{???}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{digits}{Set the digits for printing when \code{verbose = TRUE}.}
  \item{family}{The \pkg{bamlss} family object, see \code{\link{family.bamlss}}.
    In this case this is the \code{jm.bamlss} family object.}
  \item{start}{??? Unterschied zu jm.start?}
  \item{jm.start}{A named numeric vector containing possible starting values from a previous fit, e.g. by \code{\link{jm.mode}}.}
  \item{n.iter}{the number of MCMC iterations.}
  \item{burnin}{the burnin phase of the sampler, i.e., the number of starting samples that
    should be removed.}
  \item{thin}{the thinning parameter for MCMC simulation. E.g., \code{thin = 10} means,
    that only every 10th sampled parameter will be stored.}
  \item{step}{How many times should algorithm runtime information be printed, divides \code{n.iter}.}
}



\details{
We refer to the paper of Koehler, Umlauf, and Greven (2016) for details on the flexible additive joint model. In short, we model the hazard of subject \eqn{i} an event at time \eqn{t} as 
\deqn{h_{i}(t)= \exp [\eta_{\lambda i}(t)+
\eta_{\gamma i}+\eta_{\alpha i}(t) \cdot \eta_{\mu i}(t)]}
with predictor \eqn{\eta_{\lambda}} for all survival covariates that are time-varying or have a time-varying coefficient (including the log baseline hazard), predictor \eqn{\eta_{\gamma}} for baseline survival covariates, predictor \eqn{\eta_{\alpha}} representing the potentially time-varying association between the longitudinal marker \eqn{\eta_{\mu}} and the hazard. The longitudinal response \eqn{y_{ij}} at time points \eqn{t_{ij}} is modelled as \deqn{y_{ij}=\eta_{\mu i}(t_{ij})+e_{ij}} with independent normal errors \eqn{N(0, \exp[\eta_{\sigma i}(t_{ij})]^2)}. 

Each predictor \eqn{\eta_{ki}} is a structured additive predictor, i.e. a sum of functions of covariates \eqn{\eta_{ki} =  \sum_{m=1}^{M_k} f_{km}(\bm{x}_{ki})}. Each of these functions can be modelled parametrically or using basis function evaluations from the smooth constructors in \pkg{mgcv} such as \code{\link[mgcv]{s}}, 
    \code{\link[mgcv]{te}} and \code{\link[mgcv]{ti}} and can include smooth time-varying, random or spatial effects. For the Bayesian estimation of these effects we specify corresponding priors: For linear or parametric terms we use vague normal priors, smooth and random effect terms are regularized by placing generic multivariate normal priors on the coefficients and for anisotropic smooths, when multiple smoothing variance parameters are involved, more complex prior are in place (cf. Koehler, Umlauf & Greven, 2016). We use inverse Gamma hyperpriors, i.e. IG(0.001, 0.001) to obtain an inverse Gamma full conditional for the variance parameters. We estimate the posterior mode by maximizing the log-posterior of the model using a Newton-Raphson procedure, the posterior mean is obtained via derivative-based Metropolis-Hastings sampling. We recommend to use posterior mode estimates for a quick model assessment. In order to draw correct inferences from the model, posterior mean estimates should be computed. 
We approximate integration in the survival part of the likelihood using trapezoidal rule. For posterior mode estimation 
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
# with random intercepts/slopes and a linear effect of time,
# constant association alpha and no effect of the derivative
d <- simJM(nsub = 200, long_setting = "linear", alpha_setting = "constant",
            dalpha_setting = "zero", full = FALSE)

## Formula of the according joint model
f <- list(Surv2(survtime, event, obs = y) ~ s(survtime, bs = "ps"),
          gamma ~ s(x1, bs = "ps"),
          mu ~ obstime + s(id, bs = "re")+ s(id, obstime, bs = "re"),
          sigma ~ 1,
          alpha ~ 1,
          dalpha ~ -1)

## Joint model estimation
# The family object jm.bamlss sets the default optimizer and sampler function.
# First, posterior mode estimates are computed using function jm.mode(), 
# afterwards the sampler jm.mcmc() is started.
b <- bamlss(f, data = d, family = "jm", timevar = "obstime", idvar = "id")

## Plot estimated effects.
plot(b)



## Simulate survival data
# with functional random intercepts and a nonlinear effect of time,
# time-varying association alpha and no effect of the derivative.
# Note: This specification is the default in the data simulation.
d <- simJM(nsub = 200, full = FALSE)

## Formula of the according joint model
# specify the number of knots for the smooth nonlinear effect of time
klong <- 8
f <- list(Surv2(survtime, event, obs = y) ~ s(survtime, bs = "ps"),
          gamma ~ s(x1, bs = "ps"),
          mu ~ ti(id, bs = "re") + ti(obstime, bs = "ps", k = klong) + 
               ti(id, obstime, bs = c("re", "ps"), k = c(nlevels(d$id), klong)) + 
               s(x2, bs = "ps"),
          sigma ~ 1,
          alpha ~ s(survtime, bs = "ps"),
          dalpha ~ -1)

## Estimating posterior mode only using jm.mode()
b_mode <- bamlss(f, data = d, family = "jm", timevar = "obstime", idvar = "id", 
                 sampler = MVNORM)

## Estimating posterior means using jm.mcmc() and using mode starting values
b_mean <- bamlss(f, data = d, family = "jm", timevar = "obstime", idvar = "id",
                 optimizer = FALSE, jm.start = parameters(b_mode))
plot(b_mean, model = "alpha")
plot(b_mean, model = "alpha", which = samples)
}

\keyword{regression}

