### R code from vignette source 'R2BayesX.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width = 70, prompt = "R> ", continue = "+  ")
set.seed(1090)
library("R2BayesX")
data("ZambiaBnd")
data("BeechBnd")
options("use.akima" = TRUE)


###################################################
### code chunk number 2: data-illustration
###################################################
data("ZambiaNutrition", "ZambiaBnd", package = "R2BayesX")


###################################################
### code chunk number 3: fit-illustration (eval = FALSE)
###################################################
## b <- bayesx(stunting ~ sx(agechild) + sx(mbmi) +
##   sx(district, bs = "gk", map = ZambiaBnd),
##   family = "gaussian", method = "MCMC", data = ZambiaNutrition)


###################################################
### code chunk number 4: cache-illustration
###################################################
if(file.exists("illustration-model.rda")) {
load("illustration-model.rda")
} else {
b <- bayesx(stunting ~ sx(agechild) + sx(mbmi) +
  sx(district, bs = "gk", map = ZambiaBnd),
  family = "gaussian", method = "MCMC", data = ZambiaNutrition)
save(b, file = "illustration-model.rda")
}


###################################################
### code chunk number 5: summary-illustration
###################################################
summary(b)


###################################################
### code chunk number 6: plot-illustration-mbmi
###################################################
par(mar = c(4.1, 4.1, 0.1, 1.1))
plot(b, term = "sx(mbmi)")


###################################################
### code chunk number 7: plot-illustration-agechild
###################################################
par(mar = c(4.1, 4.1, 0.1, 1.1))
plot(b, term = "sx(agechild)", residuals = TRUE, cex = 0.1, rug = FALSE)


###################################################
### code chunk number 8: plot-illustration-district
###################################################
par(mar = c(0, 0, 0, 0))
plot(b, term = "sx(district)", map = ZambiaBnd, swap = TRUE, pos = "topleft")


###################################################
### code chunk number 9: illustration-plot-mbmi (eval = FALSE)
###################################################
## plot(b, term = "sx(mbmi)")


###################################################
### code chunk number 10: illustration-plot-agechild (eval = FALSE)
###################################################
## plot(b, term = "sx(agechild)", residuals = TRUE)


###################################################
### code chunk number 11: summary-illustration (eval = FALSE)
###################################################
## plot(b, term = "sx(district)", map = ZambiaBnd)


###################################################
### code chunk number 12: implementation-bayesx.construct
###################################################
bayesx.construct(sx(x, bs = "ps"))


###################################################
### code chunk number 13: bayesx.term.options1 (eval = FALSE)
###################################################
## bayesx.term.options(bs = "ps", method = "MCMC")


###################################################
### code chunk number 14: bayesx.term.options2
###################################################
out <- capture.output(bayesx.term.options(bs = "ps", method = "MCMC"))
writeLines(c(out[1:9], "..."))


###################################################
### code chunk number 15: data-zambia
###################################################
data("ZambiaNutrition", package = "R2BayesX")


###################################################
### code chunk number 16: data-zambia-bnd
###################################################
data("ZambiaBnd", package = "R2BayesX")


###################################################
### code chunk number 17: plot-zambia-map-01 (eval = FALSE)
###################################################
## plot(ZambiaBnd)


###################################################
### code chunk number 18: plot-zambia-map-02
###################################################
par(mar = c(0, 0, 0, 0))
plot(ZambiaBnd, col = "lightgray")


###################################################
### code chunk number 19: formula-zambia
###################################################
f <- stunting ~ memployment + urban + gender + meducation + sx(mbmi) +
  sx(agechild) + sx(district, bs = "mrf", map = ZambiaBnd) +
  sx(district, bs = "re")


###################################################
### code chunk number 20: fit-zambia-model (eval = FALSE)
###################################################
## zm <- bayesx(f, family = "gaussian", method = "MCMC", iterations = 12000,
##   burnin = 2000, step = 10, seed = 123, data = ZambiaNutrition)


###################################################
### code chunk number 21: cache-zambia-model
###################################################
if(file.exists("zambia-model.rda")) {
load("zambia-model.rda")
} else {
zm <- bayesx(f, family = "gaussian", method = "MCMC", iterations = 12000,
  burnin = 2000, step = 10, seed = 123, data = ZambiaNutrition)
save(zm, file = "zambia-model.rda")
}


###################################################
### code chunk number 22: summary-zambia-model
###################################################
summary(zm)


###################################################
### code chunk number 23: zambia-agechild-mbmi-plot (eval = FALSE)
###################################################
## plot(zm, term = c("sx(mbmi)", "sx(agechild)"))


###################################################
### code chunk number 24: zambia-mbmi
###################################################
par(mar = c(4.1, 4.1, 0.1, 1.1))
plot(zm, term = "sx(mbmi)")


###################################################
### code chunk number 25: zambia-agechild
###################################################
par(mar = c(4.1, 4.1, 0.1, 1.1))
plot(zm, term = "sx(agechild)")


###################################################
### code chunk number 26: zambia-district-example-kde (eval = FALSE)
###################################################
## plot(zm, term = c("sx(district):mrf", "sx(district):re"))


###################################################
### code chunk number 27: zambia-district-structured-kde
###################################################
par(mar = c(4.1, 4.1, 0.4, 1.1))
plot(zm, term = "sx(district):mrf", map = FALSE, main = "")


###################################################
### code chunk number 28: zambia-district-unstructured-kde
###################################################
par(mar = c(4.1, 4.1, 0.4, 1.1))
plot(zm, term = "sx(district):re", map = FALSE, main = "")


###################################################
### code chunk number 29: zambia-district-example (eval = FALSE)
###################################################
## plot(zm, term = "sx(district):mrf", map = ZambiaBnd)


###################################################
### code chunk number 30: zambia-district-structured
###################################################
par(mar = c(0, 0, 0, 0))
plot(zm, term = "sx(district):mrf", map = ZambiaBnd, swap = TRUE, pos = "topleft")


###################################################
### code chunk number 31: zambia-district-unstructured-samescale
###################################################
par(mar = c(0, 0, 0, 0))
plot(zm, term = "sx(district):re", map = ZambiaBnd, swap = TRUE,
  range = c(-0.32, 0.32), lrange = c(-0.32, 0.32), pos = "topleft")


###################################################
### code chunk number 32: zambia-district-example-redraw (eval = FALSE)
###################################################
## plot(zm, term = "sx(district):re", map = ZambiaBnd,
##   range = c(-0.32, 0.32), lrange = c(-0.32, 0.32))


###################################################
### code chunk number 33: zambia-mbmi-coef-samples-do
###################################################
par(oma = c(0.01, 0.01, 0.01, 0.01))
plot(zm, term = "sx(mbmi)", which = "coef-samples", main = NA)


###################################################
### code chunk number 34: zambia-autocorr-01 (eval = FALSE)
###################################################
## plot(zm, term = "sx(mbmi)", which = "var-samples", acf = TRUE)


###################################################
### code chunk number 35: zambia-mbmi-coef-samples (eval = FALSE)
###################################################
## plot(zm, term = "sx(mbmi)", which = "coef-samples")


###################################################
### code chunk number 36: zambia-autocorr-02 (eval = FALSE)
###################################################
## plot(zm, which = "max-acf")


###################################################
### code chunk number 37: zambia-autocorr-03
###################################################
par(mar = c(4.1, 4.1, 0.1, 1.1))
plot(zm, term = "sx(mbmi)", which = "var-samples", acf = TRUE, main = "")


###################################################
### code chunk number 38: zambia-autocorr-04
###################################################
par(mar = c(4.1, 4.1, 0.1, 1.1))
plot(zm, which = "max-acf", main = "")


###################################################
### code chunk number 39: fit-zambia-model-2chains (eval = FALSE)
###################################################
## zm2 <- bayesx(f, family = "gaussian", method = "MCMC", iterations = 12000,
##   burnin = 2000, step = 10, seed = 123, data = ZambiaNutrition,
##   chains = 2)


###################################################
### code chunk number 40: cache-zambia-model
###################################################
if(file.exists("zambia-model-2chains.rda")) {
load("zambia-model-2chains.rda")
} else {
zm2 <- bayesx(f, family = "gaussian", method = "MCMC", iterations = 12000,
  burnin = 2000, step = 10, seed = 123, data = ZambiaNutrition,
  chains = 2)
save(zm2, file = "zambia-model-2chains.rda")
}


###################################################
### code chunk number 41: zambia-samples2
###################################################
zs <- samples(zm2, term = "linear-samples")


###################################################
### code chunk number 42: zambia-samples2
###################################################
gelman.diag(zs, multivariate = TRUE)


###################################################
### code chunk number 43: forest-model-formula-01
###################################################
data("ForestHealth", package = "R2BayesX")
f <- defoliation ~  stand + fertilized + humus + moisture + alkali + ph +
  soil + sx(age) + sx(inclination) + sx(canopy) + sx(year) + sx(elevation)


###################################################
### code chunk number 44: fit-forest-model-01 (eval = FALSE)
###################################################
## fm1 <- bayesx(f, family = "cumlogit", method = "REML",
##   data = ForestHealth)


###################################################
### code chunk number 45: fit-forest-model-02 (eval = FALSE)
###################################################
## data("BeechBnd", package = "R2BayesX")
## fm2 <- update(fm1, . ~ . +
##   sx(id, bs = "gs", map = BeechBnd, nrknots = 20))


###################################################
### code chunk number 46: cache-forest-model
###################################################
if(file.exists("forest-model.rda")) {
load("forest-model.rda")
} else {
fm1 <- bayesx(f, family = "cumlogit", method = "REML",
  data = ForestHealth)
data("BeechBnd", package = "R2BayesX")
fm2 <- update(fm1, . ~ . +
  sx(id, bs = "gs", map = BeechBnd, nrknots = 20))
save(fm1, fm2, file = "forest-model.rda")
}


###################################################
### code chunk number 47: fit-forest-model-01-plots (eval = FALSE)
###################################################
## plot(fm1, term = c("sx(age)", "sx(inclination)", "sx(canopy)", "sx(year)",
##   "sx(elevation)"))


###################################################
### code chunk number 48: forest-no-spatial-age
###################################################
par(mar = c(4.1, 4.1, 0.1, 2.1))
plot(fm1, term = "sx(age)")


###################################################
### code chunk number 49: forest-no-spatial-inclination
###################################################
par(mar = c(4.1, 4.1, 0.1, 2.1))
plot(fm1, term = "sx(inclination)")


###################################################
### code chunk number 50: forest-no-spatial-canopy
###################################################
par(mar = c(4.1, 4.1, 0.1, 2.1))
plot(fm1, term = "sx(canopy)")


###################################################
### code chunk number 51: forest-no-spatial-year
###################################################
par(mar = c(4.1, 4.1, 0.1, 2.1))
plot(fm1, term = "sx(year)")


###################################################
### code chunk number 52: forest-no-spatial-elevation
###################################################
par(mar = c(4.1, 4.1, 0.1, 2.1))
plot(fm1, term = "sx(elevation)")


###################################################
### code chunk number 53: fit-forest-model-02-show (eval = FALSE)
###################################################
## data("BeechBnd", package = "R2BayesX")
## fm2 <- update(fm1, . ~ . +
##   sx(id, bs = "gs", map = BeechBnd, nrknots = 20))


###################################################
### code chunk number 54: summary-forest-model
###################################################
BIC(fm1, fm2)
GCV(fm1, fm2)


###################################################
### code chunk number 55: summary-forest-model
###################################################
summary(fm1)
summary(fm2)


###################################################
### code chunk number 56: forest-spatial-inclination
###################################################
par(mar = c(4.1, 4.1, 0.1, 2.1))
plot(fm2, term = "sx(inclination)")


###################################################
### code chunk number 57: forest-spatial-elevation
###################################################
par(mar = c(4.1, 4.1, 0.1, 2.1))
plot(fm2, term = "sx(elevation)")


###################################################
### code chunk number 58: forest-spatial-age
###################################################
par(mar = c(4.1, 4.1, 0.1, 2.1))
plot(fm2, term = "sx(age)")


###################################################
### code chunk number 59: forest-spatial-id (eval = FALSE)
###################################################
## plot(fm2, term = "sx(id)", map = FALSE)


###################################################
### code chunk number 60: forest-spatial-id-restrict-kde
###################################################
par(mar = c(4.1, 4.1, 0.4, 1.1))
plot(fm2, term = "sx(id)", map = FALSE, main = "")


###################################################
### code chunk number 61: forest-spatial-id
###################################################
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(fm2, term = "sx(id)", map = BeechBnd,
  height = 0.07, width = 0.27, pos = "topleft")


###################################################
### code chunk number 62: forest-spatial-id-restrict
###################################################
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(fm2, term = "sx(id)", map = BeechBnd,
  height = 0.07, width = 0.27,
  interp = TRUE, outside = TRUE,
  p.cex = 0.46, pos = "topleft")


###################################################
### code chunk number 63: forest-spatial-id (eval = FALSE)
###################################################
## plot(fm2, term = "sx(id)", map = BeechBnd)


###################################################
### code chunk number 64: forest-spatial-id-restrict (eval = FALSE)
###################################################
## plot(fm2, term = "sx(id)", map = BeechBnd,
##   interp = TRUE, outside = TRUE)


###################################################
### code chunk number 65: fit-zambia-model-step-01 (eval = FALSE)
###################################################
## f <- stunting ~ memployment + urban + gender +
##   sx(meducation, bs = "factor") + sx(mbmi) + sx(agechild) +
##   sx(district, bs = "mrf", map = ZambiaBnd) + sx(district, bs = "re")
## zms <- bayesx(f, family = "gaussian", method = "STEP",
##   algorithm = "cdescent1", startmodel = "empty", seed = 123,
##   data = ZambiaNutrition)


###################################################
### code chunk number 66: fit-zambia-model-step-02 (eval = FALSE)
###################################################
## zmsccb <- bayesx(f, family = "gaussian", method = "STEP",
##   algorithm = "cdescent1", startmodel = "empty", CI = "MCMCselect",
##   iterations = 10000, step = 10, seed = 123, data = ZambiaNutrition)


###################################################
### code chunk number 67: fit-zambia-model-step-03 (eval = FALSE)
###################################################
## zmsccb2 <- bayesx(f, family = "gaussian", method = "STEP",
##   CI = "MCMCbootstrap", bootstrapsamples = 99, iterations = 10000,
##   step = 10, seed = 123, data = ZambiaNutrition)


###################################################
### code chunk number 68: cache-zambia-model-step
###################################################
if(file.exists("zambia-model-step.rda")) {
load("zambia-model-step.rda")
} else {
data("ZambiaNutrition", "ZambiaBnd", package = "R2BayesX")
f <- stunting ~ memployment + urban + gender +
  sx(meducation, bs = "factor") + sx(mbmi) + sx(agechild) +
  sx(district, bs = "mrf", map = ZambiaBnd) + sx(district, bs = "re")
zms <- bayesx(f, family = "gaussian", method = "STEP",
  algorithm = "cdescent1", startmodel = "empty", seed = 123,
  data = ZambiaNutrition)
zmsccb <- bayesx(f, family = "gaussian", method = "STEP",
  algorithm = "cdescent1", startmodel = "empty", CI = "MCMCselect",
  iterations = 10000, step = 10, seed = 123, data = ZambiaNutrition)
zmsccb2 <- bayesx(f, family = "gaussian", method = "STEP",
  CI = "MCMCbootstrap", bootstrapsamples = 99, iterations = 10000,
  step = 10, seed = 123, data = ZambiaNutrition)
save(zms, zmsccb, zmsccb2, file = "zambia-model-step.rda")
}


###################################################
### code chunk number 69: zambia-model-step-summary
###################################################
summary(zms)


###################################################
### code chunk number 70: fit-zambia-model-step-02-show (eval = FALSE)
###################################################
## zmsccb <- bayesx(f, family = "gaussian", method = "STEP",
##   algorithm = "cdescent1", startmodel = "empty", CI = "MCMCselect",
##   iterations = 10000, step = 10, seed = 123, data = ZambiaNutrition)


###################################################
### code chunk number 71: zambia-model-step-summary-2
###################################################
summary(zmsccb)


###################################################
### code chunk number 72: fit-zambia-model-step-03-show (eval = FALSE)
###################################################
## zmsccb2 <- bayesx(f, family = "gaussian", method = "STEP",
##   CI = "MCMCbootstrap", bootstrapsamples = 99, iterations = 10000,
##   step = 10, seed = 123, data = ZambiaNutrition)


###################################################
### code chunk number 73: zambia-model-step-freqs
###################################################
term.freqs(zmsccb2, term = "sx(mbmi)")


###################################################
### code chunk number 74: fit-zambia-model-step-05 (eval = FALSE)
###################################################
## f <- stunting ~ memployment + urban + gender +
##   sx(meducation, bs = "factor") + sx(mbmi, dfstart = 2) +
##   sx(district, bs = "mrf", map = ZambiaBnd, dfstart = 5) +
##   sx(district, bs = "re", dfstart = 5) + sx(agechild, dfstart = 2)


###################################################
### code chunk number 75: fit-zambia-model-step-06 (eval = FALSE)
###################################################
## zmsud <- bayesx(f, family = "gaussian", method = "STEP",
##   algorithm = "cdescent1", startmodel = "userdefined", CI = "MCMCselect",
##   iterations = 10000, step = 10, seed = 123, data = ZambiaNutrition)


###################################################
### code chunk number 76: appendix-s1
###################################################
bayesx.construct(sx(mbmi))


###################################################
### code chunk number 77: appendix-s2
###################################################
bayesx.construct(s(mbmi, bs = "ps"))


###################################################
### code chunk number 78: appendix-te (eval = FALSE)
###################################################
## sx(mbmi, agechild, bs = "te")
## te(mbmi, agechild, bs = "ps", k = 7)


###################################################
### code chunk number 79: appendix-mrf (eval = FALSE)
###################################################
## sx(district, bs = "mrf", map = ZambiaBnd)
##  s(district, bs = "mrf", xt = list(map = ZambiaBnd))


###################################################
### code chunk number 80: large-data (eval = FALSE)
###################################################
## set.seed(321)
## file <- paste(tempdir(), "/data.raw", sep = "")
## n <- 5e+06
## dat <- data.frame(x = rep(runif(1000, -3, 3), length.out = n))
## dat$y <- with(dat, sin(x) + rnorm(n, sd = 2))
## write.table(dat, file = file, quote = FALSE, row.names = FALSE)


###################################################
### code chunk number 81: large-data-01 (eval = FALSE)
###################################################
## b <- bayesx(y ~ sx(x), family = "gaussian", method = "MCMC",
##   iterations = 3000, burnin = 1000, step = 2, predict = FALSE,
##   data = file, seed = 123)


