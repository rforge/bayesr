\name{jm.mcmc}
\alias{jm.mcmc}

\title{
  Joint model posterior mean estimation using MCMC
}

\description{
  This sampler function implements a derivative based Markov Chain Mone Carlo (MCMC) algorithm for
  flexible joint models with structured additive predictors.
}

\usage{
jm.mcmc(x, y, family, start, jm.start = NULL, weights, offset,
  n.iter = 1200, burnin = 200, thin = 1,
  verbose = TRUE, digits = 4, step = 20,
  fixed = NULL, slice = NULL, ...)
}

\arguments{
 \item{x}{The \code{x} list, as returned from function
    \code{\link{bamlss.frame}} and transformed by function \code{\link{jm.transform}},
    holding all model matrices and other information that is used for
    fitting the model.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{family}{A \pkg{bamlss} family object, see \code{\link{family.bamlss}}.
    In this case this is the \code{\link{cox.bamlss}} family object.}
  \item{start}{A named numeric vector containing possible starting values, the names are based on
    function \code{\link{parameters}}. Korrekt?}
  \item{jm.start}{A named numeric vector containing possible starting values from a previous fit, e.g. by \code{\link{jm.mode}}.}
  \item{weights}{Prior weights on the data, as returned from function \code{\link{bamlss.frame}}.}
  \item{offset}{Can be used to supply model offsets for use in fitting,
    returned from function \code{\link{bamlss.frame}}.}
  \item{n.iter}{Sets the number of MCMC iterations.}
  \item{burnin}{Sets the burnin phase of the sampler, i.e., the number of starting samples that
    should be removed.}
  \item{thin}{Defines the thinning parameter for MCMC simulation. E.g., \code{thin = 10} means,
    that only every 10th sampled parameter will be stored.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{digits}{Set the digits for printing when \code{verbose = TRUE}.}
  \item{step}{How many times should algorithm runtime information be printed, divides \code{n.iter}.}
  \item{\dots}{Currently not used.}
}

\details{
  The sampler uses derivative based proposal functions to create samples of parameters.
  The proposals are based on one Newton-Raphson iteration centered
  at the last state. The integrals that are part of the time-dependent function updates are solved numerically. In addition, smoothing variances are sampled using slice sampling.
}

\value{
  The function returns samples of parameters. The samples are provided as a
  \code{\link[coda]{mcmc}} matrix.
}

\references{
  Umlauf N, Klein N, Zeileis A (2016). Bayesian Additive Models for Location
    Scale and Shape (and Beyond). \emph{(to appear)}
}

\seealso{
  \code{\link{cox.mcmc}}, \code{\link{cox.bamlss}}, \code{\link{surv.transform}},
  \code{\link{simSurv}}, \code{\link{bamlss}}
}

\examples{
\dontrun{
set.seed(123)

## Simulate survival data
# with random intercepts and linear effect of time,
# constant association alpha and no effect of the derivative
d <- sim_jm(###)

## Formula of the according joint model
f <- list(Surv2(survtime, event, obs = y) ~ s(survtime, bs="ps"),
          gamma ~ s(x1, bs="ps"),
          mu ~ obstime + s("id", bs="re")
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
# time-varying association alpha and no effect of the derivative
d <- sim_jm(###)

## Formula of the according joint model
# specify degrees of freedom for the smooth nonlinear effect of time
long_df <- 6
f <- list(Surv2(survtime, event, obs = y) ~ s(survtime, bs="ps"),
          gamma ~ s(x1, bs="ps"),
          mu ~ ti(id, bs="re") + ti(obstime, bs="ps", k=nlevels(d$id)) + 
               ti(id, obstime, bs=c("re", "ps"), k=c(nlevels(d$id), long_df)) + s(x2, bs="ps")
          sigma ~ 1,
          alpha ~ s(survtime, bs="ps"),
          dalpha ~ -1)

## Estimating with jm.mode() only
b_mode <- bamlss(f, data=d, family="jm", timevar="obstime", idvar="id", 
             sampler = MVNORM)

## Estimating with jm.mcmc() and using mode starting values
b_mean <- bamlss(f, data=d, family="jm", timevar="obstime", idvar="id",
                 optimizer=FALSE, jm.start = parameters(b_mode))
}
}

\keyword{regression}
\keyword{survival}

