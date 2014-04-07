### R code from vignette source 'R2BayesX.Rnw'
set.seed(1090)
library("R2BayesX")
data("ZambiaNutrition", "ZambiaBnd", package = "R2BayesX")

### Section 2 - Motivating example
## model formula of the illustration model
f <- stunting ~ sx(agechild) + sx(mbmi) +
  sx(district, bs = "gk", map = ZambiaBnd, full = TRUE)

## fitting the model with bayesx()
b <- bayesx(f, family = "gaussian", method = "MCMC",
  data = ZambiaNutrition)

## summary statistics
summary(b)

## plot of the effect of the mother's bmi
plot(b, term = "sx(mbmi)")

## plot of the effect of the age of the child
## with partial residuals
plot(b, term = "sx(agechild)", residuals = TRUE, cex = 0.1, rug = FALSE)

## plot of the spatial effect using a map of Zambia
plot(b, term = "sx(district)", map = ZambiaBnd, swap = TRUE)


### Section 4.2 - Model specification
## corresponding BayesX command using the constructor
## function bayesx.construct()
bayesx.construct(sx(x, bs = "ps"))


### Section 5.2 - Available additive terms
## show possible options for various types of bases and methods  
bayesx.term.options(bs = "ps", method = "MCMC")


### Section 6.1 - Childhood malnutrition in Zambia: Analysis with MCMC
## plot of the map of Zambia
plot(ZambiaBnd, col = "lightgray")

## the formula of the Zambia model
f <- stunting ~ memployment + urban + gender + meducation + sx(mbmi) +
  sx(agechild) + sx(district, bs = "mrf", map = ZambiaBnd) +
  sx(district, bs = "re")

## estimate the model using MCMC
zm <- bayesx(f, family = "gaussian", method = "MCMC", iterations = 12000,
  burnin = 2000, step = 10, seed = 123, data = ZambiaNutrition)

## summary statistics
summary(zm)

## plot of the effect of the mother's bmi
plot(zm, term = "sx(mbmi)")

## plot of the effect of the age of the child
plot(zm, term = "sx(agechild)")

## kernel density estimate of the structured spatial effect
plot(zm, term = "sx(district):mrf", map = FALSE)

## kernel density estimate of the unstructured spatial effect
plot(zm, term = "sx(district):re", map = FALSE)

## map effect plot of the structured spatial effect
plot(zm, term = "sx(district):mrf", map = ZambiaBnd, swap = TRUE)

## map effect plot of the unstructured spatial effect
## using the same range for the colors as for the
## structured spatial effect
plot(zm, term = "sx(district):re", map = ZambiaBnd, swap = TRUE,
  range = c(-0.32, 0.32), lrange = c(-0.32, 0.32))

## sampling paths of the coefficients of the P-spline
## term for the mother's bmi 
plot(zm, term = "sx(mbmi)", which = "coef-samples", ask = TRUE)

## autocorrelation of sampling paths of the variance parameter
## of the P-spline term for the mother's bmi 
plot(zm, term = "sx(mbmi)", which = "var-samples", acf = TRUE)

## plot of the maximum autocorrelation of all parameters
## of the Zambia model
plot(zm, which = "max-acf", main = "")

## estimate model using 2 chains
zm2 <- bayesx(f, family = "gaussian", method = "MCMC", iterations = 12000,
  burnin = 2000, step = 10, seed = 123, data = ZambiaNutrition, chains = 2)

## extract parametric parameter samples 
zs <- samples(zm2, term = "linear-samples")

## compute the Gelman's and Rubin convergence diagnostic
gelman.diag(zs, multivariate = TRUE)


### Section 6.2 - Forest health dataset: Analysis with REML
data("ForestHealth", "BeechBnd", package = "R2BayesX")

## the model formula of the example
f <- defoliation ~  stand + fertilized + humus + moisture + alkali + ph +
  soil + sx(age) + sx(inclination) + sx(canopy) + sx(year) + sx(elevation)

## first model, without a spatial effect
fm1 <- bayesx(f, family = "cumlogit", method = "REML",
  data = ForestHealth)

## second model including a spatial effect
fm2 <- update(fm1, . ~ . +
  sx(id, bs = "gs", map = BeechBnd, nrknots = 20))

## plot of the estimated effects of the first model
plot(fm1, term = "sx(age)")
plot(fm1, term = "sx(inclination)")
plot(fm1, term = "sx(canopy)")
plot(fm1, term = "sx(year)")
plot(fm1, term = "sx(elevation)")

## the resukting BICs and GCV scores of the two models
BIC(fm1, fm2)
GCV(fm1, fm2)

## summary statistics
summary(fm1)
summary(fm2)

## plots of the estimated effects of the second model
plot(fm2, term = "sx(inclination)")
plot(fm2, term = "sx(elevation)")
plot(fm2, term = "sx(age)")

## kernel density estimate of the spatial effect
plot(fm2, term = "sx(id)", map = FALSE, main = "")

## map effect plot of the spatial effect
plot(fm2, term = "sx(id)", map = BeechBnd,
  height = 0.24, width = 0.41)

## map effect plot of the spatial effect
plot(fm2, term = "sx(id)", map = BeechBnd,
  height = 0.07, width = 0.27, pos = "topleft")

## interpolated spatial effect plot
plot(fm2, term = "sx(id)", map = BeechBnd,
  height = 0.07, width = 0.27,
  interp = TRUE, outside = TRUE,
  p.cex = 0.46, pos = "topleft")


### Section 6.3 - Childhood malnutrition in Zambia: Analysis with STEP
## the stepwise model formula
f <- stunting ~ memployment + urban + gender +
  sx(meducation, bs = "factor") + sx(mbmi) + sx(agechild) +
  sx(district, bs = "mrf", map = ZambiaBnd) + sx(district, bs = "re")

## estimating the model using the stepwise algorithm
zms <- bayesx(f, family = "gaussian", method = "STEP",
  algorithm = "cdescent1", startmodel = "empty", seed = 123,
  data = ZambiaNutrition)

## estimating the model using the stepwise algorithm
## including confidence intervals computed via simulation
zmsccb <- bayesx(f, family = "gaussian", method = "STEP",
  algorithm = "cdescent1", startmodel = "empty", CI = "MCMCselect",
  iterations = 10000, step = 10, seed = 123, data = ZambiaNutrition)

##  summary statistics, with and without confidence intervals
summary(zms)
summary(zmsccb)

## model with unconditioal confidence bands
zmsccb2 <- bayesx(f, family = "gaussian", method = "STEP",
  CI = "MCMCbootstrap", bootstrapsamples = 99, iterations = 10000,
  step = 10, seed = 123, data = ZambiaNutrition)

## show selection frequencies
term.freqs(zmsccb2, term = "sx(mbmi)")


### Appendic C - Smooth term constructor functions
## examples of different defaults and resulting
## paramaters using functions s() of package mgcv
## and the main model term constructor function sx()
## in R2BayesX
bayesx.construct(sx(mbmi))
bayesx.construct(s(mbmi, bs = "ps"))

