## Required packages.
library("bamlss")
library("gamlss.dist")
library("gamlss.tr")

## Load the data.
data("FlashAustria", package = "FlashAustria")

## Generate zero truncated Sichel distribution.
ztSICHEL <- trun(0, family = "SICHEL", local = FALSE)

## Generate zero truncated beta negative binomial.
ztBNB <- trun(0, family = "BNB", local = FALSE)

## Model formula, up to four parameters.
f <- list(
  counts ~ s(sqrt_cape, bs = "ps") + s(hcc, bs = "ps") + s(str, bs = "ps") + s(tciw, bs = "ps") + s(q_prof_PC1, bs = "ps"),
         ~ s(sqrt_cape, bs = "ps") + s(hcc, bs = "ps") + s(str, bs = "ps") + s(tciw, bs = "ps") + s(q_prof_PC1, bs = "ps"),
         ~ s(sqrt_cape, bs = "ps") + s(hcc, bs = "ps") + s(str, bs = "ps") + s(tciw, bs = "ps") + s(q_prof_PC1, bs = "ps"),
         ~ s(sqrt_cape, bs = "ps") + s(hcc, bs = "ps") + s(str, bs = "ps") + s(tciw, bs = "ps") + s(q_prof_PC1, bs = "ps")
)

## Estimate models.
set.seed(123)
b_ztnbinom <- bamlss(f, data = FlashAustriaTrain,       ## Standard interface.
  family = "ztnbinom", binning = TRUE,                  ## General arguments.
  optimizer = opt_boost, maxit = 1000,                  ## Boosting arguments.
  thin = 50, burnin = 2000, n.iter = 3000, cores = 50,  ## Sampler arguments.
  light = TRUE)

set.seed(123)
b_ztSICHEL <- bamlss(f, data = FlashAustriaTrain,
  family = ztSICHEL, binning = TRUE,
  optimizer = opt_boost, maxit = 1000,
  thin = 50, burnin = 2000, n.iter = 3000, cores = 50,
  light = TRUE)

set.seed(123)
b_ztBNB <- bamlss(f, data = FlashAustriaTrain,
  family = ztBNB, binning = TRUE,
  optimizer = opt_boost, maxit = 1000,
  thin = 50, burnin = 2000, n.iter = 3000, cores = 50,
  light = TRUE)

save(b_ztnbinom, b_ztSICHEL, b_ztBNB, file = "FlashAustriaModel.rda")

## Final model.
set.seed(123)
b_ztSICHEL_f <- bamlss(f, data = FlashAustriaTrain,
  family = ztSICHEL, binning = TRUE,
  optimizer = opt_boost, maxit = 1000,
  thin = 200, burnin = 2000, n.iter = 6000, cores = 50,
  light = TRUE)

save(b_ztSICHEL_f, file = "FLashAustiaModel_f.rda")

