\name{r}
\alias{r}
\alias{random}

\title{
  Construct Hierarchical and Non-Hierarchical Random Effects Model Terms
}

\description{
  Function used in definition of random terms within \code{\link{bayesx}} model formulae. The 
  function does not evaluate design matrices etc. - it exists purely to help set up a model using 
  random effect terms.
}

\usage{
r(x, h = NULL, by = NA, xt = NULL, data = NULL, 
  weights = NULL, subset = NULL, offset = NULL, 
  na.action = na.fail, contrasts = NULL, 
  control = bayesx.control(...), ...)
}
 
\arguments{
  \item{x}{the factor/cluster or group variable for which random effects should be estimated. If
    \code{method = "HMCMC"} in \code{\link{bayesx}}, \code{x} may also be a 
    \code{\link[stats]{formula}} if hierarchical covariates are available that are used for the 
    model, please see argument \code{h} and the examples.}
  \item{h}{optional, only if \code{method = "HMCMC"} in \code{\link{bayesx}}, argument 
    \code{h} may specify a second level equation that is used for hierarchical model 
    estimation. I.e., in a hierarchical setting, certain variables are only measured on a cluster 
    or group level, hence, different observations are only given by the number of different clusters 
    or groups in the data. Utilizing this hierarchical structure will increase sampling 
    time and mixing. Suppose the following data structure: a group variable \eqn{id} with \eqn{n} 
    observations indicating in which of the \eqn{k} groups observation \eqn{i} lies in. Secondly: 
    two variables \eqn{x_1} and \eqn{x_2}, that are only measured within the group level, i.e. the 
    number of different observations for \eqn{x_1} and \eqn{x_2} is \eqn{k}. Also suppose these 
    variables may have a smooth and a linear effect on the response \eqn{y} 
    (\eqn{y = (y_1,\ldots,y_n)}). Now, the hierarchical model is specified with 
    \code{y ~ r(id, method = ~ s(x1) + x2)}. Note that within the second equation, 
    \code{~ s(x1) + x2}, all model terms specified in \code{\link{bayesx.construct}} may be applied, 
    i.e. also random effects terms with higher levels may enter the second level equation. Please 
    see the examples for further illustrations.}
  \item{by}{specify a variable for which the slope estimator is modeled as random.}
  \item{xt}{specify extra information for estimation.}
  \item{data, weights, subset, offset, na.action, contrasts, control, \dots}{for hierarchical models
    , additional model specifications may be supplied here, please see the examples below, also see 
    \code{\link{bayesx}}.}
}

\value{
  A \code{list} is returned where information on labels and variables and other, depending on the 
  specification, is provided. 
}

\note{
  For hierarchical models, plotting of the total effect is initialized by an extra argument
  \code{total = TRUE}, please see the last example.
}

\author{
  Thomas Kneib, Stefan Lang, Nikolaus Umlauf, Achim Zeileis.
}

\seealso{
  \code{\link{bayesx}}, \code{\link{plot.bayesx}}.
}

\examples{
r(id1 ~ x1 + s(x2) + r(id2 ~ s(x3) + x4))

## same with
r(id1, ~ x1 + s(x2) + r(id2 ~ s(x3) + x4))

\dontrun{
## example for use within bayesx()
## estimate simple random effects model
set.seed(111)
n <- 1000
     
## index vector for random effects
N <- 100
dat <- data.frame(id = sort(rep(1:N, n/N)), 
  x1 = runif(n, -3, 3))
     
## create some iid normal random 
##  effects 
dat$re <- with(dat, rnorm(N, sd = 0.6)[id])

## response
dat$y <- with(dat, 1.5 + sin(x1) + re + rnorm(n, sd = 0.6))
     
## estimate model
b <- bayesx(y ~ s(x1, bs = "ps") + r(id), data = dat)
summary(b)

## extract effects and plot 
## against true random effects
## together with effect
## of x1
re.fit <- fitted(b, term = "r(id)")[,2L]
op <- par(no.readonly = TRUE)
par(mfrow=c(1,2))
ylim <- range(c(re.fit, dat$re))
plot(re.fit ~ unique(dat$re), ylim = ylim, xlim = ylim,
     xlab = "True effects", ylab = "Fitted random effects")
abline(a = 0, b = 1)
plot(b, term = "s(x1)")
par(op)


## now estimate a random effects
## model with a stage 2 covariate
## effect
set.seed(333)
n <- 1000

## 1st stage
N <- 100
dat1 <- data.frame(id = sort(rep(1:N, n/N)), 
  x1 = runif(n, -3, 3))
dat1$re <- with(dat1, rnorm(N, sd = 0.6)[id])

## 2nd stage
dat2 <- data.frame(id = unique(dat1$id), 
  x2 = runif(N, -3, 3))

## response
dat1$y <- with(dat1, 1.5 + sin(x1) + re + cos(dat2$x2)[id] + rnorm(n, sd = 0.6))

## estimate hierarchical model
## with the intercept in the 
## 2nd stage
system.time(b1 <- bayesx(y ~ -1 + s(x1, bs = "ps") + 
  r(id ~ 1 + s(x2 , bs = "ps"), data = dat2), 
  method = "HMCMC", data = dat1))
summary(b1)

## non hierarchical
dat1$x2 <- with(dat1, dat2$x2[id])
system.time(b2 <- bayesx(y ~ s(x1, bs = "ps") + s(x2, bs = "ps") + r(id),
  method = "MCMC", data = dat1))
summary(b2)

## plot x2 effects
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
plot(b1, term = "s(x2)", main = "Hierarchical")
plot(b2, term = "s(x2)", main = "Non Hierarchical")
par(op)


## now a spatial example
## with structured and
## unstructered spatial 
## effect, the structured
## effect is modeled in level 2
set.seed(333)

## simulate some geographical data
data("MunichBnd")
N <- length(MunichBnd); names(MunichBnd) <- 1:N
n <- N*5

## regressors
dat1 <- data.frame(id = sort(rep(1:N, n/N)), x1 = runif(n, -3, 3))
dat2 <- data.frame(re = rnorm(N, sd = 0.6), 
  sp = sort(runif(N, -2, 2), decreasing = TRUE),
  id = unique(dat1$id))

## response
dat1$y <- with(dat1, 1.5 + sin(x1) + dat2$sp[id] + dat2$re[id] + rnorm(n, sd = 0.6))

## estimate models with
## bayesx MCMC and REML
b <- bayesx(y ~ -1 + s(x1, bs = "ps") + 
  r(id ~ 1 + s(id, bs = "mrf", xt = list(map = MunichBnd)), data = dat2), 
  method = "HMCMC", data = dat1)

## plot the spatial effects
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 3))
plot(b, term = "s(id)", map = MunichBnd, 
  main = "Structured spatial effect")
plot(b, term = "r(id)", map = MunichBnd, 
  main = "Unstructured spatial effect")
plot(b, term = "r(id)", map = MunichBnd,
  main = "Total spatial effect", total = TRUE)
par(op)
}
}

\keyword{regression} 
