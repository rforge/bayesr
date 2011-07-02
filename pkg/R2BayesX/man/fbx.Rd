\name{fbx}
\alias{fbx}

\title{
  Construct BayesX Model Terms in A Formula
}

\description{
  Function \code{fbx} is a model term constructor function for terms used within the \code{formula}
  argument of function \code{\link{bayesx}}. The function does not evaluate matrices etc., the 
  behavior is similar to function \code{\link[mgcv]{s}} from package \code{\link[mgcv]{mgcv}}. It 
  purely exists to build a basic setup for the model term which can be processed by function 
  \code{\link{bayesx.construct}}.
}

\usage{
fbx(x, z = NA, bs = "ps", by = NA, ...)
}
 
\arguments{
  \item{x}{the covariate the model term is a function of.}
  \item{z}{a second covariate.}
  \item{bs}{a \code{\link{character}} string, specifying the basis/type which is used for this model 
    term.}
  \item{by}{a \code{\link{numeric}} or \code{\link{factor}} variable of the same dimension as each 
    covariate. In the \code{numeric} vector case the elements multiply the smooth, evaluated at the 
    corresponding covariate values (a `varying coefficient model' results). In the \code{factor} 
    case the term is replicated for each \code{factor} level. Note that centering of the term may be 
    needed, please see the notes.}
  \item{\dots}{special controlling arguments or objects used for the model term, see also 
    the examples and function \code{\link{bayesx.term.options}} for all possible optional 
    parameters.}
}

\details{
  The following term types may be specified using argument \code{bs}:
  \itemize{
    \item \code{"lasso"}, \code{"nigmix"}, \code{"ridge"}: shrinkage of fixed effects, defines a 
      shrinkage-prior for the corresponding parameters \eqn{\gamma_j}, \eqn{j = 1, \ldots, q}, 
      \eqn{q \geq 1} of the linear effects \eqn{x_1, \ldots, x_q}. There are three priors possible: 
      ridge-, lasso- and Normal Mixture of inverse Gamma (NMIG)-prior. 
    \item \code{"bl"}, \code{"baseline"}: nonlinear baseline effect in hazard regression or 
      multi-state models, defines a P-spline with second order random walk penalty for the 
      parameters of the spline for the log-baseline effect \eqn{log(\lambda(time))}.
    \item \code{"kr"}, \code{kriging}: kriging with stationary Gaussian random fields. Covariate
      arguments \code{x} and \code{z} must be specified.
    \item \code{"gk"} or \code{"geokriging"}: geokriging with stationary Gaussian random fields: 
      estimates a stationary Gaussian random field based on the centroids of a map object provided 
      in boundary format, see function \code{\link[BayesX]{read.bnd}}, within the \code{xt} argument 
      of function \code{\link[mgcv]{s}}, e.g. \code{xt = list(map = MapBnd)}. 
    \item \code{"gs"}, \code{"geospline"}: geosplines based on two-dimensional P-splines with 
      first order random walk penalty, defines a two-dimensional P-spline for the spatial covariate 
      region with a two-dimensional first order random walk penalty for the parameters of the 
      spline. Estimation is based on the coordinates of the centroids of the regions provided by a 
      map object in in boundary format, see function \code{\link[BayesX]{read.bnd}}, within the 
      \code{xt} argument of function \code{\link[mgcv]{s}}, e.g. \code{xt = list(map = MapBnd)}, or 
      within argument \code{"..."} of function \code{f}, i.e. \code{f(id, bs = "gs", map = MapBnd)}. 
    \item \code{"mrf"}, \code{"spatial"}: Markov random fields: defines a Markov random field 
      prior for a spatial covariate, where geographical informations are provided by a map object in
      boundary or graph file format, see function \code{\link[BayesX]{read.bnd}} and 
      \code{\link[BayesX]{read.gra}}, within the \code{xt} argument of function 
      \code{\link[mgcv]{s}}, e.g. \code{xt = list(map = MapBndorGra)}, or within argument 
      \code{"..."} of function \code{f}, i.e. \code{f(id, bs = "mrf", map = MapBnd)}.
    \item \code{"ps"}, \code{"psplinerw1"}, \code{"psplinerw2"}, \code{"pspline"}: P-spline with 
      first or second order difference penalty. 
    \item \code{"rw1"}, \code{"rw2"}: zero degree P-splines, defines a zero degree P-spline with 
      first or second order difference penalty. A zero degree P-spline typically estimates for every 
      distinct covariate value in the data set a separate parameter. Usually there is no reason to 
      prefer zero degree P-splines over higher order P-splines. An exception are ordinal covariates 
      or continuous covariates with only a small number of different values. For ordinal covariates 
      higher order P-splines are not meaningful while zero degree P-splines might be an alternative 
      to modeling nonlinear relationships via a dummy approach with completely unrestricted 
      regression parameters.
    \item \code{"season"}: seasonal effect of a time scale. 
    \item \code{"te"}, \code{"pspline2dimrw2"}: defines a two-dimensional P-spline based on the 
      tensor product of one-dimensional P-splines with a two-dimensional first order random walk 
      penalty for the parameters of the spline. Covariate arguments \code{x} and \code{z} must be 
      specified.
    \item \code{"offset"}: an offset term.
    \item \code{"ra"}, \code{"random"}: i.i.d. Gaussian (random) effect of a group indicator 
      covariate, e.g. the covariate may be an individual indicator when analyzing longitudinal data.
    \item \code{"factor"}: special BayesX specifier for \code{\link[base]{factor}}s, especially 
      meaningful if \code{method = "STEP"} in \code{\link{bayesx}}, since the factor term is then 
      treated as a full term, which is either included or removed from the model.
    \item \code{"cs"}, \code{"catspecific"}: category-specific linear effect (in cumulative or 
      sequential models only).
  } 
}

\note{
  Some care has to be taken with the identifiability of varying coefficients terms. The standard in
  BayesX is to center nonlinear main effects terms around zero whereas varying coefficient terms are
  not centered. This makes sense since main effects nonlinear terms are not identifiable and varying
  coefficients terms are usually identifiable. However, there are situations where a varying
  coefficients term is not identifiable. Then the term must be centered. Since centering is not
  automatically accomplished it has to be enforced by the user by adding option
  \code{center = TRUE} in function \code{f}. To give an example, the varying coefficient terms in 
  \eqn{\eta = \ldots + g_1(z_1)z + g_2(z_2)z + \gamma_0 + \gamma_1 z + \ldots} are not
  identified, whereas in \eqn{\eta = \ldots + g_1(z_1)z + \gamma_0 + \ldots}, the varying 
  coefficient term is identifiable. In the first case, centering is necessary, in the second case, 
  it is not.
}

\value{
  A \code{list} of class \code{"xx.smooth.spec"}, where \code{"xx"} is a basis/type identifying code
  given by the \code{bs} argument of \code{f}. 
}

\author{
  Thomas Kneib, Stefan Lang, Nikolaus Umlauf, Achim Zeileis.
}

\seealso{
  \code{\link{bayesx}}, \code{\link{bayesx.term.options}}, \code{\link[mgcv]{s}}, \code{\link{r}},
  \code{\link{bayesx.construct}}.
}

\examples{
## funktion fbx() returns a list
## which is then processed by function
## bayesx.construct to build the
## BayesX model term structure
fbx(x)

bayesx.construct(fbx(x))
bayesx.construct(fbx(x, bs = "rw1"))
bayesx.construct(fbx(x, bs = "factor"))
bayesx.construct(fbx(x, bs = "offset"))
bayesx.construct(fbx(x, z, bs = "te"))

## using a map for markov random fields
data("FantasyBnd")
plot(FantasyBnd)
bayesx.construct(fbx(id, bs = "mrf", map = FantasyBnd))

## examples using optional controlling
## parameters and objects
bayesx.construct(fbx(x, bs = "ps", knots = 20))
bayesx.construct(fbx(x, bs = "ps", nrknots = 20))
bayesx.construct(fbx(x, bs = "ps", knots = 20, nocenter = TRUE))

## use of bs with original 
## BayesX syntax 
bayesx.construct(fbx(x, bs = "psplinerw1"))
bayesx.construct(fbx(x, bs = "psplinerw2"))
bayesx.construct(fbx(x, z, bs = "pspline2dimrw2"))
bayesx.construct(fbx(id, bs = "spatial", map = FantasyBnd))
bayesx.construct(fbx(x, z, bs = "kriging"))
bayesx.construct(fbx(id, bs = "geospline", map = FantasyBnd, nrknots = 5))
bayesx.construct(fbx(x, bs = "catspecific"))


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
b1 <- bayesx(y ~ fbx(x), method = "REML", data = dat)

## note that the shortcut f() instead of fbx()
## is working within bayesx() 
b1 <- bayesx(y ~ f(x), method = "REML", data = dat)

## increase inner knots
## decrease degree of the P-spline
b2 <- bayesx(y ~ f(x, knots = 30, degree = 2), method = "REML", data = dat)


## compare reported output
summary(c(b1, b2))

## plot the effect for both models
plot(c(b1, b2), residuals = TRUE)


## more examples
set.seed(111)
n <- 500

## regressors
dat <- data.frame(x = runif(n, -3, 3), z = runif(n, -3, 3),
  w = runif(n, 0, 6), fac = factor(rep(1:10, n/10)))

## response
dat$y <- with(dat, 1.5 + sin(x) + cos(z) * sin(w) +
  c(2.67, 5, 6, 3, 4, 2, 6, 7, 9, 7.5)[fac] + rnorm(n, sd = 0.6))

## estimate model
b <- bayesx(y ~ f(x) + f(z, w, bs = "te") + fac,
  data = dat, method = "MCMC")

summary(b)
plot(b)


## now a mrf example
## note: the regional identification
## covariate and the map regionnames
## should be coded as integer
set.seed(333)
     
## simulate some geographical data
data("MunichBnd")
N <- length(MunichBnd); n <- N*5
names(MunichBnd) <- 1:N
     
## regressors
dat <- data.frame(x1 = runif(n, -3, 3),
  id = as.factor(rep(names(MunichBnd), length.out = n)))
dat$sp <- with(dat, sort(runif(N, -2, 2), decreasing = TRUE)[id])
     
## response
dat$y <- with(dat, 1.5 + sin(x1) + sp + rnorm(n, sd = 1.2))

## estimate models with
## bayesx MCMC and REML
b <- bayesx(y ~ f(x1) + f(id, bs = "mrf", map = MunichBnd), 
  method = "REML", data = dat)

## summary statistics
summary(b)

## plot the effects
op <- par(no.readonly = TRUE)
par(mfrow = c(1,2))
plot(b, term = "f(id)", map = MunichBnd, 
  main = "bayesx() estimate")
plotmap(MunichBnd, x = dat$sp, id = dat$id, 
  main = "Truth")
par(op)


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
b <- bayesx(y ~ f(x1, bs = "ps") + f(id, bs = "ra"), data = dat)
summary(b)
plot(b)

## extract estimated random effects
## and compare with true effects
plot(fitted(b, term = "f(id)")$Mean ~ unique(dat$re))
}
}

\keyword{regression} 
