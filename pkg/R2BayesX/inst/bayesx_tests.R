## comprehensive BayesX testing
library("R2BayesX")

## load the artificial data set
data("GAMart")


## MCMC
## numeric response
b <- bayesx(num ~ fac + sx(x1) + sx(x2) + sx(x3) +
  sx(long, lat, bs = "te") + sx(id, bs = "re"),
  data = GAMart, method = "MCMC", family = "gaussian")
summary(b)
plot(b)

## gamma
GAMart$num2 <- GAMart$num + abs(min(GAMart$num)) + 1

## working
b <- bayesx(num2 ~ fac + sx(x1) + sx(x2) + sx(x3) +
  sx(long, lat, bs = "te"),
  data = GAMart, method = "MCMC", family = "gamma")

## FIXME: not working?
b <- bayesx(num2 ~ fac + sx(x1) + sx(x2) + sx(x3) +
  sx(long, lat, bs = "te") + sx(id, bs = "re"),
  data = GAMart, method = "MCMC", family = "gamma")

## binomial
## FIXME: factor and random effects not working?
b <- bayesx(bin ~ fac + sx(x1) + sx(x2) + sx(x3) +
  sx(long, lat, bs = "te") + sx(id, bs = "re"),
  data = GAMart, method = "MCMC", family = "binomial")

## working
b <- bayesx(bin ~ sx(x1) + sx(x2) + sx(x3) + sx(long, lat, bs = "te"),
  data = GAMart, method = "MCMC", family = "binomial")

## FIXME: not working?
b <- bayesx(bin ~ fac + sx(x1) + sx(x2) + sx(x3) + sx(long, lat, bs = "te"),
  data = GAMart, method = "MCMC", family = "binomial")

## FIXME: not working?
b <- bayesx(bin ~ sx(x1) + sx(x2) + sx(x3) + sx(long, lat, bs = "te") + sx(id, bs = "re"),
  data = GAMart, method = "MCMC", family = "binomial")

## FIXME: not working?
b <- bayesx(bin ~ sx(x1) + sx(id, bs = "re"),
  data = GAMart, method = "MCMC", family = "binomial")

## check with mgcv gam()
b <- gam(bin ~ fac + s(x1) + s(x2) + s(x3) +
  s(long, lat) + s(id, bs = "re"),
  data = GAMart, family = binomial)
summary(b)
plot(b)

## cumprobit
## FIXME: factor and random effects not working?
b <- bayesx(cat ~ fac + sx(x1) + sx(x2) + sx(x3) +
  sx(long, lat, bs = "te") + sx(id, bs = "re"),
  data = GAMart, method = "MCMC", family = "cumprobit")

## working, but category specific thresholds are missing?
b <- bayesx(cat ~ fac + sx(x1) + sx(x2),
  data = GAMart, method = "MCMC", family = "cumprobit")
summary(b)
plot(b)

## FIXME: not working?
b <- bayesx(cat ~ fac + sx(x1) + sx(x2) + sx(x3),
  data = GAMart, method = "MCMC", family = "cumprobit")


## REML
## numeric response
b <- bayesx(num ~ fac + sx(x1) + sx(x2) + sx(x3) +
  sx(long, lat, bs = "te") + sx(id, bs = "re"),
  data = GAMart, method = "REML", family = "gaussian")
summary(b)
plot(b)

## binomial
b <- bayesx(bin ~ fac + sx(x1) + sx(x2) + sx(x3) +
  sx(long, lat, bs = "te") + sx(id, bs = "re"),
  data = GAMart, method = "REML", family = "binomial")
summary(b)
plot(b)

## cumprobit
b <- bayesx(cat ~ fac + sx(x1) + sx(x2) + sx(x3) +
  sx(long, lat, bs = "te") + sx(id, bs = "re"),
  data = GAMart, method = "REML", family = "cumprobit", outfile = "~/tmp")
summary(b)
plot(b)

