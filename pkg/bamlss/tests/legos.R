library("bamlss")
data("ZambiaNutrition", package = "R2BayesX")

## (1) Model formulae.
bamlss.formula(stunting ~ s(mbmi) + meducation)

f <- list(
  stunting ~ s(mbmi) + meducation,
  sigma ~ s(mbmi)
)

bamlss.formula(f)
bamlss.formula(f, family = gamma.bamlss())

## Automatic filling with intercept.
bamlss.formula(stunting ~ mbmi, family = zinb.bamlss())
bamlss.formula(f, family = zinb.bamlss())


## (2) Parsing model input.
pm <- parse.input.bamlss(f, data = ZambiaNutrition, family = gaussian)
names(pm)

## Note that parse.input.bamlss() may handle special user defined smooths
## in addition to mgcv user defined smooths, one just needs to add a specials = TRUE
## within an smooth.construct() call. This is usefull to e.g. estimate NURBS with JAGS.


## (3) Different engines in use.
##     These steps are all combined in bamlss()/xreg().
##     To ease formulae specifictaion within bamlss(), one only needs
##     to write e.g. zinb instead of zinb.bamlss or zinb.bamlss().

## BayesX, needs sx().
f2 <- list(
  stunting ~ sx(mbmi) + meducation,
  sigma ~ sx(mbmi)
)
pm2 <- parse.input.bamlss(f2, data = ZambiaNutrition, family = gaussian)

## Transformer to parse sx() smooths correctly.
tpm <- transformBayesX(pm2)

## Setup the model code and write out data.
sm <- setupBayesX(tpm)

## Start sampling.
ms <- samplerBayesX(sm)
plot(ms)

## Get results, returns object of class "bamlss".
mr <- resultsBayesX(tpm, ms)
plot(mr)
summary(mr)

## Now with JAGS.
tpm <- transformJAGS(pm)
sm <- setupJAGS(tpm)
ms <- samplerJAGS(sm)
mr <- resultsJAGS(tpm, ms)
plot(mr)
summary(mr)

## Try STAN, uses same transformer and results as JAGS.
## WARNING: takes very very long!
tpm <- transformJAGS(pm)
sm <- jags2stan(tpm)
ms <- samplerSTAN(sm, n.iter = 1200, burnin = 200, thin = 1)
mr <- resultsJAGS(tpm, ms)
plot(mr)
summary(mr)

## With IWLS, does not need a setup function (work in progress).
tpm <- transformIWLS(pm)
ms <- samplerIWLS(tpm, method = c("backfitting", "MCMC"), maxit = 100)
mr <- resultsIWLS(tpm, ms)
plot(mr)
summary(mr)

