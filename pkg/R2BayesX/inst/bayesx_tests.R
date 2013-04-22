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

## binomial
## FIXME: factor and random effects not working?
b <- bayesx(bin ~ fac + sx(x1) + sx(x2) + sx(x3) +
  sx(long, lat, bs = "te") + sx(id, bs = "re"),
  data = GAMart, method = "MCMC", family = "binomial")
summary(b)
plot(b)

## cumprobit
## FIXME: factor and random effects not working?
b <- bayesx(cat ~ fac + sx(x1) + sx(x2) + sx(x3) +
  sx(long, lat, bs = "te") + sx(id, bs = "re"),
  data = GAMart, method = "MCMC", family = "cumprobit")
summary(b)
plot(b)


## REML
## numeric response
b <- bayesx(num ~ fac + sx(x1) + sx(x2) + sx(x3) +
  sx(long, lat, bs = "te") + sx(id, bs = "re"),
  data = GAMart, method = "REML", family = "gaussian")
summary(b)
plot(b)

## binomial
## FIXME: factor and random effects not working?
b <- bayesx(bin ~ fac + sx(x1) + sx(x2) + sx(x3) +
  sx(long, lat, bs = "te") + sx(id, bs = "re"),
  data = GAMart, method = "REML", family = "binomial")
summary(b)
plot(b)

## cumprobit
## FIXME: factor and random effects not working?
b <- bayesx(cat ~ fac + sx(x1) + sx(x2) + sx(x3) +
  sx(long, lat, bs = "te") + sx(id, bs = "re"),
  data = GAMart, method = "REML", family = "cumprobit")
summary(b)
plot(b)
