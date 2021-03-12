\name{opt_bbfit}
\alias{bbfit}
\alias{opt_bbfit}
\alias{bbfitp}
\alias{opt_bbfitp}
\alias{contribplot}

\title{Batchwise Backfitting}

\description{
  Batchwise backfitting estimation engine for GAMLSS using very large data sets.
}

\usage{
## Batchwise backfitting engine.
opt_bbfit(x, y, family, shuffle = TRUE, start = NULL, offset = NULL,
  epochs = 1, nbatch = 10, verbose = TRUE, ...)

bbfit(x, y, family, shuffle = TRUE, start = NULL, offset = NULL,
  epochs = 1, nbatch = 10, verbose = TRUE, ...)

## Parallel version.
opt_bbfitp(x, y, family, mc.cores = 1, ...)

## Loglik contribution plot.
contribplot(x, ...)
}

\arguments{
   \item{x}{For function \code{bfit()} the \code{x} list, as returned from function
    \code{\link{bamlss.frame}}, holding all model matrices and other information that is used for
    fitting the model. For the updating functions an object as returned from function
    \code{\link[mgcv]{smooth.construct}} or \code{\link[mgcv]{smoothCon}}. For function
    \code{contribplot()}, a \code{"bamlss"} object using \code{bbfit()} with argument
    \code{select = TRUE}.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{family}{A \pkg{bamlss} family object, see \code{\link{family.bamlss}}.}
  \item{shuffle}{Should observations be shuffled?}
  \item{start}{A named numeric vector containing possible starting values, the names are based on
    function \code{\link{parameters}}.}
  \item{offset}{Can be used to supply model offsets for use in fitting,
    returned from function \code{\link{bamlss.frame}}.}
  \item{epochs}{For how many epochs should the algorithm run?}
  \item{nbatch}{Number of batches. Can also be a number between 0 and 1, i.e., determining
    the fraction of observations that should be used for fitting.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{mc.cores}{On how many cores should estimation be started?}
  \item{\dots}{For \code{bbfitp()} all arguments to be passed to \code{bbfit()}.}
}

\details{
  The algorithm uses batch-wise estimation of regression coefficients and smoothing variances.
  The smoothing variances are estimated on an hold-out batch. This way, models for very large
  data sets can be estimated. Note, the algorithm can work in combination with the \pkg{ff}
  and \pkg{ffbase} package, i.e., the entire data is never in the computer RAM. Therefore, the
  data can either to be stored as comma separated file on disc or provided as \code{"ffdf"}
  data frame, see also the examples.

  The optimizer functions use additional arguments:
  \itemize{
    \item \code{batch_ids}. This argument can either be a list of indices specifying the
      batches that should be used for estimation, or a vector of length 2, where the
      first element specifies the number of observations that should be sampled for each
      batch and the second argument specifies the number of batches, see the example.
    \item \code{nu}, the step length control parameter. Defaults to \code{nu = 0.05}.
      If argument \code{slice = TRUE} then \code{nu = 1}.
    \item \code{loglik}, defaults to \code{loglik = FALSE}. If set to \code{loglik = TRUE}
      the "out-of-sample" log-likelihood is used for smoothing variance estimation.
    \item \code{aic}, defaults to \code{aic = FALSE}, If set to \code{aic = TRUE}
      the "out-of-sample" AIC is used for smoothing variance estimation.
    \item \code{eps_loglik}, defaults to \code{eps_loglik = 0.01}. This argument specifies
      the relative change in the "out-of-sample" log-likelihood that is needed such that
      a model term gets updated.
    \item \code{select}, defaults to \code{select = FALSE}. If set to \code{select = TRUE},
      the algorithm only selects the model term with the largest contribution in the
      "out-of-sample" log-likelihood for updating in each iteration/batch.
    \item \code{always}, defaults to \code{always = FALSE}. If set to \code{always = TRUE}
      no log-likelihood contribution checks will be used and model terms are always updated.
    \item \code{K}, defaults to \code{K = 2}. This argument controls the penalty on the
      degrees of freedom in the computation of the AIC.
    \item \code{slice}, defaults to \code{slice = FALSE}. If set to \code{slice = TRUE},
      slice sampling using the "out-of-sample" log-likelihood or AIC is used for smoothing
      variance estimation. Moreover, \code{always = TRUE}, \code{eps_loglik = -Inf} and
      \code{nu = 1}.
  }

  When using function \code{opt_bbfitp}, the parameter updates are stored as \code{"mcmc"}
  objects. In this case the traceplots can be visualized using \code{\link{plot.bamlss}}.
}

\value{
  For function \code{opt_bbfit()} a list containing the following objects:
  \item{fitted.values}{A named list of the fitted values of the modeled parameters
    of the selected distribution.}
  \item{parameters}{The estimated set regression coefficients and smoothing variances.}
  \item{shuffle}{Logical}
  \item{runtime}{The runtime of the algorithm.}
}

\seealso{
\code{\link{bamlss}}, \code{\link{bfit}}
}

\examples{
\dontrun{## Simulate data.
set.seed(123)
d <- GAMart(n = 27000, sd = -1)

## Write data to disc.
tf <- tempdir()
write.table(d, file.path(tf, "d.raw"), quote = FALSE, row.names = FALSE, sep = ",")

## Model formula.
f <- list(
  num ~ s(x1,k=40) + s(x2,k=40) + s(x3,k=40) + te(lon,lat,k=10),
  sigma ~ s(x1,k=40) + s(x2,k=40) + s(x3,k=40) + te(lon,lat,k=10)
)

## Different flavors:
## (1) Using out-of-sample aic for smoothing
##     variance estimation. Update is only accepted
##     if the "out-of-sample" log-likelihood is
##     increased.

## Specify batches.
batch_ids <- c("nobs" = 2000, "nbatch" = 20)

b1 <- bamlss(f, data = d, sampler = FALSE, optimizer = opt_bbfit,
  batch_ids = batch_ids, nu = 0.1, aic = TRUE, eps_loglik = -Inf,
  always = FALSE)

plot(b1)

## (2) Same but always update, this mimics the classic SGD.
##     Note that we use opt_bbfitp() to inspect the coefficient
##     paths. Note, for prediction only the last iterations should
##     be used.
b2 <- bamlss(f, data = d, sampler = FALSE, optimizer = opt_bbfitp,
  batch_ids = batch_ids, nu = 0.1, aic = TRUE, eps_loglik = -Inf,
  always = TRUE)

plot(b2, which = "samples")
}
}

\keyword{regression}

