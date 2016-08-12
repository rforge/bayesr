\name{jm.transform}
\alias{jm.transform}

\title{
  Joint Model Transformer Function
}

\description{
  This function takes a \code{\link{bamlss.frame}} and computes design matrices of model terms
  based on a time grid for time-dependent structured additive predictors in a survival context.
  Note that this transformer function is usually used internally by function \code{\link{bamlss}}
  and is the default transformer function using the \code{\link{jm.bamlss}} family object.

  The time grid design matrices can be used to construct the full structured additive predictor
  for each time point. This way it is possible to solve the integrals that are part of, e.g.,
  a Newton-Raphson updating scheme, numerically.

  See the example section on how to extract the time grid design matrices.
}

\usage{
jm.transform(x, y, data, terms, knots, formula, family, jm.start = NULL,
  subdivisions = 25, timedependent = c("lambda", "mu", "alpha", "dalpha"), 
  timevar = NULL, idvar = NULL, alpha = .Machine$double.eps, mu = NULL, sigma = NULL,
  sparse = TRUE, ...)
}

\arguments{
 \item{x}{The \code{x} list, as returned from function \code{\link{bamlss.frame}},
    holding all model matrices and other information that is used for
    fitting the model.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{data}{The data.frame that should be used for setting up all matrices.}
  \item{terms}{List of terms (one list element for each predictor) from the formula setup}
  \item{knots}{Potential specification of knots placed in one of the terms. Structure???}
  \item{formula}{List of formulas (one list element for each predictor) as returned from function
    \code{\link{bamlss.frame}}}
  \item{family}{A \pkg{bamlss} family object, see \code{\link{family.bamlss}}.
    In this case this is the \code{\link{jm.bamlss}} family object.}
  \item{jm.start}{Starting parameters for the estimation. If \code{NULL} ???. See Details on how to set starting parameters}
  \item{subdivisions}{Number of time points created for each individual for numerical integration.}
  \item{timedependent}{A character vector specifying the names of parameters in \code{x}
    that are time-dependent. Time grid design matrices are only computed for these parameters.}
  \item{timevar}{A character specifying the name of the longitudinal time variable in the
    data set.}
  \item{idvar}{Name of the variable specifying the identifier of individuals. 
  Identifier needs to be a factor.}
  \item{alpha}{A value for the intercept of a parameter names alpha. Typically the association
    parameter of a longitudinal and survival process in a joint model.???}
  \item{mu}{???}
  \item{sigma}{Wird nicht benutzt. Raus?}
  \item{sparse}{Logical. Should sparse setup be used setting up design matrices for \code{mu}}
  \item{\dots}{Arguments passed to function \code{\link{bamlss.engine.setup}}.}
}

\value{
  A \code{\link{bamlss.frame}} including the time grid design matrices.
}

\details{
The function performs various steps to set up the design matrices for a 
joint model fit. Among these are: 
(1) create a subject-specific timegrid of integration points, 
(2) set up design matrices for the time-varying predictors for integration, 
(3) set up sparse matrices for terms within predictor \code{mu}, 
(4) remove intercept from predictor \code{lambda} (global intercept is estimated in \code{gamma}),
and (5) set up \code{eta_timegrid}, the full linear combination of all time-varying  terms, evaluated
at the timegrid.
}

\seealso{
  \code{\link{jm.bamlss}}, \code{\link{jm.mode}}, \code{\link{jm.mcmc}},
  \code{\link{jm.sim}}, \code{\link{bamlss}}
}

\examples{
set.seed(111)

## Simulate survival data
# with random intercepts and linear effect of time,
# constant association alpha and no effect of the derivative
d <- sim_jm(###)

## Formula of the joint model
f <- list(Surv2(survtime, event, obs = y) ~ s(survtime, bs="ps"),
          gamma ~ s(x1, bs="ps"),
          mu ~ obstime + s("id", bs="re"),
          sigma ~ 1,
          alpha ~ 1,
          dalpha ~ -1)

## Create the bamlss.frame
bf <- bamlss.frame(f, data=d, family="jm", timevar="obstime", idvar="id")
print(bf)

## Apply the transformer ??? (throws error)
bf <- with(bf, jm.transform(x = x, y = y, terms = terms, knots = knots,
           formula = formula, family = family, data = model.frame,
           timevar = "obstime", idvar = "id"))

}

\keyword{regression}
\keyword{survival}

