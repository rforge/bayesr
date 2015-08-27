library("bamlss")
data("GAMart", package = "R2BayesX")

## (1) Model formulae.
bamlss.formula(num ~ s(x1) + id)

f <- list(
  num ~ s(x1) + s(x3) + id,
  sigma ~ s(x2),
  id ~ s(x3)
)

bamlss.formula(f)
bamlss.formula(f, family = gamma.bamlss())

## Automatic filling with intercept.
bamlss.formula(num ~ x1, family = zinb.bamlss())
bamlss.formula(f, family = zinb.bamlss())

## Terms object from bamlss.formula().
terms(bamlss.formula(f))
terms.bamlss(f)
terms(bamlss.formula(f), model = c(1, 1), sterms = FALSE)
terms(bamlss.formula(f), model = c(1, 1), sterms = FALSE, pterms = FALSE)


## (2) Get the bamlss frame.
bf <- bamlss.frame(f, data = GAMart, family = "gaussian")
names(bf)
print(bf)
head(model.frame(bf))
bf$model.frame <- NULL
head(model.frame(bf))
model.response(model.frame(bf))
response.name(bf)
response.name(model.frame(bf))


## (3) model.matrix() and smooth.construct()
model.matrix(bf)
head(model.matrix(bf, model = c(1, 1)))
smooth.construct(bf)
bf$model.frame <- NULL
head(model.matrix(bf, model = c(1, 1)))
str(smooth.construct(bf))
str(smooth.construct(bf, model = c(1, 1)))

bf <- bamlss.frame(f, data = GAMart, family = "gaussian",
  model.matrix = FALSE, smooth.construct = FALSE)
print(bf)
model.matrix(bf)
head(model.matrix(bf, model = c(1, 1)))
smooth.construct(bf)
names(smooth.construct(bf, model = c(1, 1)))


## (4) Complex multilevel structures.
f <- list(
  cat ~ s(x1) + id,
  id ~ s(x3)
)

bf <- bamlss.frame(f, data = GAMart, family = "multinomial")
print(bf)
head(model.response(model.frame(bf)))

terms(bf)
terms(formula(bf))
terms(bf, model = c(2, 1))
terms(bf, model = c(2, 1), sterms = FALSE)
terms(bf, model = c(2, 1), pterms = FALSE)

d <- GAMart[1:30, ]

model.matrix(formula(bf), data = d)
head(model.matrix(formula(bf), data = d, model = c(1, 1)))
head(model.matrix(formula(bf, model = c(1, 1)), data = d))
head(model.matrix(terms(bf, model = c(1, 1), drop = FALSE), data = d))

smooth.construct(formula(bf), data = d)
smooth.construct(terms(bf), data = d)
smooth.construct(formula(bf), data = d, model = c(1, 1))
smooth.construct(formula(bf, model = c(1, 1)), data = d)
smooth.construct(terms(bf, model = c(1, 1), drop = FALSE), data = d)

## Extract or initiallize parameters.
p <- parameters(bf)

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

