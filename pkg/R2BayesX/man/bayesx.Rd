\name{bayesx}
\alias{bayesx}
\alias{R2BayesX}

\title{
  Estimate STAR Models with BayesX
}

\description{
  Fit structured additive regression models via Markov chain Monte Carlo simulation techniques
  (MCMC), mixed model methodology (REML) and a penalized least squares (respectively penalized 
  likelihood) approach for estimating structured additive regression models using model selection 
  tools (STEP).
}

\usage{  
bayesx(formula, data, weights = NULL, subset = NULL, 
  offset = NULL, na.action = na.fail, contrasts = NULL, 
  control = bayesx.control(\dots), \dots)
}

\arguments{
  \item{formula}{symbolic description of the model (of type \code{y ~ x}).}
  \item{data}{a \code{\link{data.frame}} or \code{\link{list}} containing the model response 
    variable and covariates required by the formula. By default the variables are taken from 
    \code{environment(formula)}: typically the environment from which \code{bayesx} is called. 
    Argument \code{data} may also be a character string defining the directory the data is stored, 
    where the first row in the data set must contain the variable names and columns should be tab
    separated.}
  \item{weights}{prior weights on the data.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting 
    process.}
  \item{offset}{can be used to supply a model offset for use in fitting.}
  \item{na.action}{a function which indicates what should happen when the data contain \code{NA}'s.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg} of 
    \code{\link[stats]{model.matrix.default}}.}
  \item{control}{specify several global control parameters for \code{bayesx}, see 
    \code{\link{bayesx.control}}.}
  \item{\dots}{arguments passed to \code{\link{bayesx.control}}.}
}

\details{  
  Estimation of regression models with \pkg{BayesX} is based on two different inferential concepts: 
  Markov chain Monte Carlo simulation techniques corresponding to full Bayesian inference 
  and mixed model methodology corresponding to penalized likelihood or empirical Bayes inference.
  Both concepts have been implemented as follows:
  
  \bold{MCMC} simulation techniques: A fully Bayesian interpretation of structured additive 
  regression models is obtained by specifying prior distributions for all unknown parameters. 
  Estimation can be facilitated using Markov chain Monte Carlo simulation techniques, a general and 
  versatile concept for Bayesian inference. Implementation of MCMC methodology thereby 
  provides numerically efficient schemes for structured additive regression models. Suitable 
  proposal densities have been developed to obtain rapidly mixing, well-behaved sampling schemes 
  without the need for manual tuning.

  \bold{REML}, mixed model based estimation: An increasingly popular way to estimate semi-parametric 
  regression models is the representation of penalization approaches as mixed models. Within 
  \pkg{BayesX} this concept has been extended to structured additive regression models and several 
  types of non-standard regression situations. The general idea is to take advantage of the close 
  connection between penalty concepts and corresponding random effects distributions. The smoothing 
  parameters of the penalties then transform to variance components in the random effects (mixed) 
  model. While the selection of smoothing parameters has been a difficult task for a long time, 
  several estimation procedures for variance components in mixed models are already available since 
  the 1970's. The most popular one is restricted maximum likelihood in Gaussian mixed models with 
  marginal likelihood as the non-Gaussian counterpart. REML employs mixed model methodology for the 
  estimation of structured additive regression models. While regression coefficients are estimated 
  based on penalized likelihood, restricted maximum likelihood or marginal likelihood estimation 
  forms the basis for the determination of smoothing parameters. From a Bayesian perspective, this 
  yields empirical Bayes/posterior mode estimates for the structured additive regression models. 
  However, estimates can also merely be interpreted as penalized likelihood estimates from a 
  frequentist perspective.
  
  \bold{STEP}, penalized least squares including model selection: As a third alternative 
  \pkg{BayesX} provides a penalized least squares (respectively penalized likelihood) approach for 
  estimating structured additive regression models. In addition to the previously described 
  estimation alternatives, a powerful variable and model selection tool is included. Model choice 
  and estimation of the parameters is done simultaneously. The algorithms are able to
  \itemize{
    \item decide whether a particular covariate enters the model,
    \item decide whether a continuous covariate enters the model linearly or nonlinearly,
    \item decide whether a spatial effect enters the model,
    \item decide whether a unit- or cluster specific heterogeneity effect enters the model
    \item select complex interaction effects (two dimensional surfaces, varying coefficient terms)
    \item select the degree of smoothness of nonlinear covariate, spatial or cluster specific
      heterogeneity effects.
    }
  Inference is based on penalized likelihood in combination with fast algorithms for selecting 
  relevant covariates and model terms. Different models are compared via various goodness of fit 
  criteria, e.g. AIC, BIC, GCV and 5 or 10 fold cross validation.
  
  Now, to estimate STAR models, the wrapper function \code{bayesx} starts by simply translating
  \R syntax into \pkg{BayesX} commands, that have originally been implemented, with function
  \code{\link{parse.bayesx.input}} and \code{\link{write.bayesx.input}}. Afterwards the computed 
  program file is send as a whole to the command-line binary executable version of \pkg{BayesX} with 
  \code{\link{run.bayesx}}, which is processed fully externally. As a last step, function 
  \code{\link{read.bayesx.output}} will read the estimated model object returned from \pkg{BayesX} 
  back into \R. 
  
  For estimation of STAR models, function \code{bayesx} thereby uses formula syntax as provided in 
  package \code{\link[mgcv]{mgcv}} (see \code{\link[mgcv]{formula.gam}}), i.e. models including 
  smooth terms may be specified using either
  \itemize{
    \item \code{\link[mgcv]{s}} type or
    \item \code{\link[mgcv]{te}} type model terms.
    }
  In addition, another term type for random effects has been implemented, also see \code{\link{r}},
  with possibilities to estimate hierarchical STAR models. For a detailed description of the syntax 
  used within \code{bayesx} models see \code{\link{bayesx.construct}} and 
  \code{\link{bayesx.term.options}}.

  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"bayesx"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[graphics]{plot}}, 
  \code{\link[stats]{residuals}} and \code{\link[stats]{fitted}}.
  
  See \code{\link{fitted.bayesx}}, \code{\link{plot.bayesx}}, and \code{\link{summary.bayesx}} for 
  more details on these methods.
}

\value{ 
  See function \code{\link{read.bayesx.output}}.
}

\note{
  If a model is specified with a structured and an unstructured spatial effect, e.g. the model 
  formula is something like \code{y ~ s(id, bs = "mrf", map = mymap) + r(id)}, the model output 
  contains of one additional total spatial effect, named with \code{"s(id):total"}. Also see the 
  last example.
}

\author{
  Thomas Kneib, Stefan Lang, Nikolaus Umlauf, Achim Zeileis.
}

\section{WARNINGS}{
  Function \code{bayesx} or \code{\link{run.bayesx}} cannot be processed if the \code{bin} path to
  the binary executable command-line version is not specified probably, i.e. before using 
  \code{bayesx} the binary has to be installed, see function \code{\link{install.bayesx}} for
  installation details. After installation the \code{bin} path can bet set with

  \code{options(bayesx.bin = "/path/to/BayesX")}

  or within function \code{\link{bayesx.control}}, where \code{"/path/to/BayesX"} represents the 
  full path to the binary, e.g. on windows systems this may be something like 
  \code{"C:/programs/BayesX/commandline/bayesx.exe"}. 

  For geographical effects, note that BayesX may crash if the region identification covariate is
  a \code{\link{factor}}, it is recommended to code these variables as \code{\link{integer}}, please
  see the example below.
}

\references{
  For methodological and reference details see the \pkg{BayesX} manuals available at:
  
  \url{http://www.stat.uni-muenchen.de/~bayesx/}.
}

\seealso{
  \code{\link{parse.bayesx.input}}, \code{\link{write.bayesx.input}}, \code{\link{run.bayesx}}, 
  \code{\link{read.bayesx.output}}, \code{\link{summary.bayesx}}, \code{\link{plot.bayesx}},
  \code{\link{fitted.bayesx}}, \code{\link{bayesx.construct}}, \code{\link{bayesx.term.options}},
  \code{\link[mgcv]{formula.gam}}, \code{\link{r}}.
}

\examples{
\dontrun{
## generate some data
set.seed(111)
n <- 200

## regressor
dat <- data.frame(x = runif(n, -3, 3))

## response
dat$y <- with(dat, 1.5 + sin(x) + rnorm(n, sd = 0.6))

## estimate models with
## bayesx REML and MCMC
b1 <- bayesx(y ~ s(x, bs = "ps"), 
  method = "REML", data = dat)
b2 <- bayesx(y ~ s(x, bs = "ps"), 
  method = "MCMC", iter = 1200, data = dat)

## compare reported output
summary(c(b1, b2))

## plot the effect for both models
plot(c(b1, b2), residuals = TRUE)

## use confint
confint(b1, level = 0.99)
confint(b2, level = 0.99)


## more examples
set.seed(111)
n <- 500

## regressors
dat <- data.frame(x = runif(n, -3, 3), z = runif(n, -3, 3),
  w = runif(n, 0, 6), fac = factor(rep(1:10, n/10)))

## response
dat$y <- with(dat, 1.5 + sin(x) + cos(z) * sin(w) +
  c(2.67, 5, 6, 3, 4, 2, 6, 7, 9, 7.5)[fac] + rnorm(n, sd = 0.6))

## estimate models with
## bayesx MCMC and REML
## and compare with
## mgcv gam()
b1 <- bayesx(y ~ s(x, bs = "ps") + s(z, w, bs = "te") + fac,
  data = dat, method = "MCMC")
b2 <- bayesx(y ~ s(x, bs = "ps") + s(z, w, bs = "te") + fac,
  data = dat, method = "REML")
b3 <- gam(y ~ s(x, bs = "ps") + te(z, w, bs = "ps") + fac, 
  data = dat)

## summary statistics
summary(b1)
summary(b2)
summary(b3)

## plot the effects
op <- par(no.readonly = TRUE)
par(mfrow = c(3,2))
plot(b1, term = "s(x)")
plot(b1, term = "s(z,w)")
plot(b2, term = "s(x)")
plot(b2, term = "s(z,w)")
plot(b3, select = 1)
vis.gam(b3, c("z","w"), theta = 40, phi = 40)
par(op)

## combine models b1 and b2
b <- c(b1, b2)

## summary
summary(b)

## only plot effect 2 of both models
plot(b, term = "s(z,w)") 

## with residuals
plot(b, term = "s(z,w)", residuals = TRUE) 

## same model with kriging
b <- bayesx(y ~ s(x, bs = "ps") + s(z, w, bs = "kr") + fac, 
  method = "REML", data = dat)
plot(b)


## now a mrf example
## note: the regional identification
## covariate and the map regionnames
## should be coded as integers
set.seed(333)

## simulate some geographical data
data("MunichBnd")
N <- length(MunichBnd); names(MunichBnd) <- 1:N
n <- N*5

## regressors
dat <- data.frame(id = rep(1:N, n/N), x1 = runif(n, -3, 3))
dat$sp <- with(dat, sort(runif(N, -2, 2), decreasing = TRUE)[id])

## response
dat$y <- with(dat, 1.5 + sin(x1) + sp + rnorm(n, sd = 0.6))

## estimate models with
## bayesx MCMC and REML
b1 <- bayesx(y ~ s(x1, bs = "ps") + 
  s(id, bs = "mrf", xt = list(map = MunichBnd)), 
  method = "MCMC", data = dat)
b2 <- bayesx(y ~ s(x1, bs = "ps") + 
  s(id, bs = "mrf", xt = list(map = MunichBnd)), 
  method = "REML", data = dat)

## summary statistics
summary(b1)
summary(b2)

## plot the effects
op <- par(no.readonly = TRUE)
par(mfrow = c(1,3))
plot(b1, term = "s(id)", map = MunichBnd, 
  main = "bayesx() MCMC estimate")
plot(b2, term = "s(id)", map = MunichBnd, 
  main = "bayesx() REML estimate")
plotmap(MunichBnd, x = unique(cbind(dat$id, dat$sp)), 
  main = "Truth")
par(op)

## try geosplines instead
b <- bayesx(y ~ s(id, bs = "gs", k = 20, xt = list(map = MunichBnd)) + 
  s(x1, bs = "ps"), data = dat)
summary(b)
plot(b, term = "s(id)", map = MunichBnd)

## geokriging
b <- bayesx(y ~ s(id, bs = "gk", k = 20, xt = list(map = MunichBnd)) + 
  s(x1, bs = "ps"), data = dat)
summary(b)
plot(b, term = "s(id)", map = MunichBnd)

## perspective plot of the effect
plot(b, term = "s(id)")

## image and contour plot 
plot(b, term = "s(id)", 
  image = TRUE, contour = TRUE, grid = 200)


## model with random effects
set.seed(333)
N <- 30
n <- N*10

## regressors
dat <- data.frame(id = sort(rep(1:N, n/N)), x1 = runif(n, -3, 3))
dat$re <- with(dat, rnorm(N, sd = 0.6)[id])

## response
dat$y <- with(dat, 1.5 + sin(x1) + re + rnorm(n, sd = 0.6))

## estimate model
b <- bayesx(y ~ s(x1, bs = "ps") + r(id), data = dat)
summary(b)
plot(b)

## extract estimated random effects
## and compare with true effects
re.fit <- fitted(b, term = "r(id)")
plot(unique(re.fit[["r(id)"]][, "Mean"]) ~ unique(dat$re))


## now a spatial example
## with structured and
## unstructered spatial 
## effect
set.seed(333)

## simulate some geographical data
data("MunichBnd")
N <- length(MunichBnd); names(MunichBnd) <- 1:N
n <- N*5

## regressors
dat <- data.frame(id = rep(1:N, n/N), x1 = runif(n, -3, 3))
dat$sp <- with(dat, sort(runif(N, -2, 2), decreasing = TRUE)[id])
dat$re <- with(dat, rnorm(N, sd = 0.6)[id])

## response
dat$y <- with(dat, 1.5 + sin(x1) + sp + re + rnorm(n, sd = 0.6))

## estimate models with
## bayesx MCMC and REML
b <- bayesx(y ~ s(x1, bs = "ps") + 
  s(id, bs = "mrf", xt = list(map = MunichBnd)) +
  r(id), method = "MCMC", data = dat)
summary(b)

## plot all spatial effects
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 3))
plot(b, term = "s(id)", map = MunichBnd, 
  main = "Structured spatial effect")
plot(b, term = "r(id)", map = MunichBnd, 
  main = "Unstructured spatial effect")
plot(b, term = "s(id):total", map = MunichBnd, 
  main = "Total spatial effect", digits = 4)
par(mfrow = op)


## some experiments with the
## stepwise algorithm
## generate some data
set.seed(321)
n <- 1000

## regressors
dat <- data.frame(x1 = runif(n, -3, 3), x2 = runif(n),
  x3 = runif(n, 3, 6), x4 = runif(n, -3, -1))

## response
dat$y <- with(dat, 1.5 + sin(x1) + 0.6 * x2 + rnorm(n, sd = 0.6))

## estimate model with STEP
b <- bayesx(y ~ s(x1, bs = "ps") + s(x2, bs = "ps") + 
  s(x3, bs = "ps") +  s(x4, bs = "ps"), 
  method = "STEP", iter = 10000, step = 10, data = dat)
summary(b)
plot(b)


## a probit example
set.seed(111)
n <- 1000
dat <- data.frame(x <- runif(n, -3, 3))

dat$z <- with(dat, sin(x) + rnorm(n))
dat$y <- rep(0, n)
dat$y[dat$z > 0] <- 1

b <- bayesx(y ~ s(x, bs = "ps"), family = "binomialprobit", data = dat)
summary(b)
plot(b)


## estimate varying coefficient models
set.seed(333)
n <- 1000
dat <- data.frame(x = runif(n, -3, 3), id = factor(rep(1:4, n/4)))

## response
dat$y <- with(dat, 1.5 + sin(x) * c(-1, 0.2, 1, 3)[id] + rnorm(n, sd = 0.6))

# estimate model
b <- bayesx(y ~ s(x, bs = "ps", by = id), method = "REML", data = dat)
summary(b)
plot(b)


## experiment with weights 
set.seed(3)
n <- 300
dat <- data.frame(x = seq(-3, 3, length = n))
dat$y <- with(dat, 1.5 + sin(x) + 
  rnorm(n, sd = (0.1 + seq(0, 0.4, length = n))))

b1 <- bayesx(y ~ s(x, bs = "ps"), data = dat)
summary(b1)
plot(b1, resid = TRUE)

dat$ressqrt <- residuals(b1)[,1L]^2
b2 <- bayesx(ressqrt ~ x, data = dat)
plot(b2, resid = TRUE)

dat$W <- exp(fitted(b2)[,1L])
b3 <- bayesx(y ~ s(x, bs = "ps"), weights = W, data = dat)

summary(c(b1, b3))
plot(c(b1, b3), resid = TRUE)
}
}

\keyword{regression}
