\name{bamlss}
\alias{bamlss}

\title{
  Fit Bayesian Additive Models for Location Scale and Shape (and Beyond)
}

\description{
  This is the main model fitting function of the package. Function \code{bamlss()}
  is a wrapper function that parses the \code{data} and the model \code{formula}, or
  extended \code{\link{bamlss.formula}}, as well as the \code{\link{bamlss.family}}
  into a \code{\link{bamlss.frame}}. The \code{\link{bamlss.frame}} then holds all model
  matrices and information that is needed for setting up estimation engines.
  The model matrices are based on \code{\link{mgcv}} infrastructures, i.e.,
  smooth terms are constructed using \code{\link[mgcv]{smooth.construct}} and
  \code{\link[mgcv]{smoothCon}}. Therefore, all \code{\link{mgcv}} model term constructors like
  \code{\link[mgcv]{s}}, \code{\link[mgcv]{te}}, \code{\link[mgcv]{t2}} and \code{\link[mgcv]{ti}}
  can be used. Identifiability conditions are imposed using function \code{\link[mgcv]{gam.side}}.

  After the \code{\link{bamlss.frame}} is set up function \code{bamlss()} applies optimizer
  and/or sampling functions using the \code{'bamlss.frame'} object. These functions can also be
  provided by the user. See the details below on how to create new engines to be used
  with function \code{bamlss()}.

  The optimum parameters and/or the samples of the parameters are then parsed to create
  summary statistics with \code{\link{summary.bamlss}}, effect plots with
  \code{\link{plot.bamlss}}, etc., see the examples.
}

\usage{
bamlss(formula, family = "gaussian", data = NULL,
  start = NULL, knots = NULL, weights = NULL,
  subset = NULL, offset = NULL, na.action = na.omit,
  contrasts = NULL, reference = NULL, transform = NULL,
  optimizer = NULL, sampler = NULL, samplestats = NULL,
  results = NULL, cores = NULL, sleep = NULL,
  combine = TRUE, model = TRUE, x = TRUE, ...)
}

\arguments{
  \item{formula}{A formula or extended formula, i.e., the \code{formula} can be a
    \code{\link{list}} of formulas where each list entry specifies the details of one parameter
    of the modeled response distribution, see \code{\link{bamlss.formula}}. For incorporating
    smooth terms, all model term constructors implemented in \code{\link{mgcv}} such as
    \code{\link[mgcv]{s}}, \code{\link[mgcv]{te}} and \code{\link[mgcv]{ti}} can be used, amongst others.}
  \item{family}{A \code{\link{bamlss.family}} object, specifying the details of the modeled
    distribution such as the parameter names, the density function, link functions, etc.}
  \item{data}{A \code{\link{data.frame}} or \code{\link{list}} containing the model
    response variable(s) and covariates specified in the \code{formula}.
    By default the variables are taken from \code{environment(formula)}:
    typically the environment from which \code{bamlss} is called.}
  \item{start}{A named numeric vector containing starting values to be send to the \code{optimizer}
    and/or \code{sampler} function. For a possible naming convention for the parameters see
    function \code{\link{parameters}}, but this is not restrictive and engine specific.}
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
  \item{transform}{A transformer function that is applied on the \code{\link{bamlss.frame}}.}
  \item{optimizer}{An optimizer function that returns, e.g., posterior mode estimates
    of the parameters as a named numeric vector. The default optimizer function is
    \code{\link{bfit}}.}
  \item{sampler}{A sampler function that returns a matrix of samples, the columns represent the
    parameters, the rows the iterations. The returned matrix must be coerced to an object of
    class \code{"mcmc"}, see \code{\link[coda]{as.mcmc}}. The default sampler function is
    \code{\link{GMCMC}}.}
  \item{samplestats}{A function computing statistics from samples, per default function
    \code{\link{samplestats}} is used.}
  \item{results}{A function computing results from the parameters and/or samples, e.g., for
    creating effect plots, see function \code{link{results.bamlss.default}}.}
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
  The main idea of this function is to provide infrastructures that make it relatively easy
  for creating estimation engines for new problems, or write interfaces to existing software
  packages, in order to estimate Bayesian additive models for location scale and shape (and beyond).

  The steps that are performed within the function are:
  \itemize{
    \item First, the function parses the \code{data}, the \code{formula} or the extended
      \code{\link{bamlss.formula}} as well as the \code{\link{bamlss.formula}} into a model frame
      like object, a \code{\link{bamlss.frame}}. This object holds all necessary model matrices
      and information that is needed for subsequent model fitting engines. Per default,
      all package \code{\link{mgcv}} smooth term constructor functions like
      \code{\link[mgcv]{s}}, \code{\link[mgcv]{te}} and \code{\link[mgcv]{ti}} can be used (see
      also function \code{\link[mgcv]{smooth.construct}}), however,
      even user defined constructors can be included, see the manual of \code{\link{bamlss.frame}}.
    \item In a second step, the \code{\link{bamlss.frame}} can be transformed, e.g., if a mixed
      model representation of smooth terms should be generated, see function \code{\link{randomize}}.
    \item Then an optimizer function is applied, e.g., a function that finds posterior mode estimates
      of the parameters. A convention for model fitting engines is that such functions should have
      the following arguments:

      \code{optimizer(x, y, family, start, weights, offset, ...)}

      Internally, function \code{bamlss()} will send the \code{x} object, the \code{y} object,
      the \code{family} object, etc., of the created \code{\link{bamlss.frame}} object to the
      optimizer function (see the manual of \code{\link{bamlss.frame}} for more details on the
      \code{x}, \code{y} and other objects). The job of the optimizer is to return a named numeric
      vector of optimum parameters. The names of the parameters should be such that they can be
      uniquely mapped to the corresponding model matrices. See function \code{\link{parameters}}
      for more details on parameter names. The default optimizer function is \code{\link{bfit}}.
  }
}

\value{
  An object of class \code{"bamlss"}.
}

\references{
%% ~put references to the literature/web site here ~
}

\author{
  Nikolaus Umlauf, Nadja Klein, Achim Zeileis.
}

\seealso{
  \code{\link{bamlss.frame}}, \code{\link{family.bamlss}}, \code{\link{formula.bamlss}},
}

\examples{
\dontrun{
## Simulated data example.
d <- GAMart()
b <- bamlss(num ~ s(x1) + s(x2) + s(x3), data = d)
summary(b)
plot(b)
}
}

\keyword{regression}
\keyword{MCMC}

