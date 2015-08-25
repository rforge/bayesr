library("bamlss")
data("GAMart", package = "R2BayesX")

## (1) Model formulae.
bamlss.formula(num ~ s(x1) + id)

f <- list(
  num ~ s(x1) + id,
  sigma ~ s(x2),
  id ~ s(x3)
)

bamlss.formula(f)
bamlss.formula(f, family = gamma.bamlss())

## Automatic filling with intercept.
bamlss.formula(num ~ x1, family = zinb.bamlss())
bamlss.formula(f, family = zinb.bamlss())


## (2) Get the bamlss frame.
bf <- bamlss.frame(f, data = GAMart, family = "gaussian")
names(bf)
head(model.frame(bf))
model.response(model.frame(bf))
response.name(bf)
response.name(model.frame(bf))

## (3) model.matrix() and smooth.construct()
model.matrix(bf)
smooth.construct(bf)

## (4) Complex multilevel structures.
f <- list(
  cat ~ s(x1) + id,
  id ~ s(x3)
)

bf <- bamlss.frame(f, data = GAMart, family = "multinomial")
names(bf$terms)
names(bf$terms$low)
head(model.response(model.frame(bf)))

terms(bf)
terms(formula(bf))
terms(bf, model = c(2, 1))

## Note that bamlss.frame() may handle special user defined smooths
## in addition to mgcv user defined smooths, one just needs to add a specials = TRUE
## within an smooth.construct() call. This is usefull to e.g. estimate NURBS with JAGS.

## (5) Run backfitting optimizer on bamlss.frame.
data("marital.nz", package = "VGAM")

bf <- bamlss.frame(mstatus ~ s(age), data = marital.nz,
  family = "multinomial", reference = "Married/Partnered")

bf <- bfit0(bf, do.optim = FALSE)

## (6) Run MCMC.
samps <- GMCMC(bf, n.iter = 1200)
plot(samps)

## (7) Process results.
b <- results(bf, samps)

## (8) Plot and summaries.
plot(b)
summary(b)

## (9) No with JAGS.
sm <- setupJAGS(bf)
samps <- samplerJAGS(sm)


f <- list(num ~ s(x1) + s(x2) + s(x3))
bf <- bamlss.frame(f, data = GAMart)
bf <- bfit0(bf)
sm <- setupJAGS(bf)
samps <- samplerJAGS(sm)



## TODOs: predict, fitted, residuals, ..., JAGS, BayesX!

