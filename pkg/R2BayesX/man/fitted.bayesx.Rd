\name{fitted.bayesx}
\alias{fitted.bayesx}

\title{
  Extract BayesX Fitted Values
}

\description{
  Extractor function to the fitted values of the estimated model with \code{\link{bayesx}} and 
  fitted model term partial effects.
}
 
\usage{
\method{fitted}{bayesx}(object, model = NULL, term = NULL, ...)
}

\arguments{
  \item{object}{an object of class \code{"bayesx"}.}
  \item{model}{for which model the fitted values should be provided, either an integer or a 
    character, e.g. \code{model = "mcmc.model"}.}
  \item{term}{if not \code{NULL}, the function will search for the term fitted values specified 
    here, either an integer or a character, eg \code{term = "s(x)"}.}
  \item{\dots}{not used.}
}

\value{ 
  Either the fitted linear predictor or if e.g. \code{term = "s(x)"}, an object with class 
  \code{"xx.bayesx"}, where \code{"xx"} is depending of the type of the term. In principle the
  returned term object is simply a matrix containing the covariate(s) and its effects, depending on 
  the estimation method, e.g. for MCMC estimated models, mean or median fitted values and other 
  quantities are returned. Several additional informations on the term are provided in the 
  \code{\link{attributes}} of the object. For all types of terms plotting functions are provided.
}

\author{
  Thomas Kneib, Stefan Lang, Nikolaus Umlauf, Achim Zeileis.
}

\seealso{
  \code{\link{read.bayesx.output}}.
}

\examples{
\dontrun{
## generate some data
set.seed(121)
n <- 500

## regressors
dat <- data.frame(x = runif(n, -3, 3), z = runif(n, 0, 1),
  w = runif(n, 0, 3))

## generate response 
dat$y <- with(dat, 1.5 + sin(x) + z -3*w + rnorm(n, sd = 0.6))

## estimate model
b1 <- bayesx(y ~ s(x, bs = "ps") + z + w, data = dat)

## extract fitted values
fit <- fitted(b1)
hist(fit, freq = FALSE)

## now extract 1st model term
## and plot it
fx <- fitted(b1, term = "s(x)")
plot(fx)
}

## now another example with
## use of read.bayesx.output
## load example data from
## package R2BayesX
dir <- paste(searchpaths()[grep("R2BayesX", searchpaths())],
             "/examples/ex02", sep = "")
b2 <- read.bayesx.output(dir)

## extract fitted values
fit <- fitted(b2)
hist(fit)

## extract model term of x1
## and plot it
fx <- fitted(b2, term = "s(x1)")
plot(fx)


## unlist and plot again
fx <- fx[["s(x1)"]]
fx[1:3,]
plot(fx)

## have a look at the attributes
names(attributes(fx))

## extract the sampling path of the variance
spv <- attr(fx, "variance.sample")
plot(spv, type = "l")


\dontrun{
## combine model objects
b <- c(b1, b2)

## extract fitted terms for second model
fit <- fitted(b, model = 2, term = 1:3)
names(fit)
plot(fit[["s(district)"]])
}
}

\keyword{regression}
