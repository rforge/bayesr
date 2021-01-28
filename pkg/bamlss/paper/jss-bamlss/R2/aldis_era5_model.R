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
  counts ~ s(sqrt_cape, bs = "ps") + s(hcc, bs = "ps") + s(str, bs = "ps") +
           s(tciw, bs = "ps") + s(q_prof_PC1, bs = "ps"),
         ~ s(sqrt_cape, bs = "ps") + s(hcc, bs = "ps") + s(str, bs = "ps") +
           s(tciw, bs = "ps") + s(q_prof_PC1, bs = "ps"),
         ~ s(sqrt_cape, bs = "ps") + s(hcc, bs = "ps") + s(str, bs = "ps") +
           s(tciw, bs = "ps") + s(q_prof_PC1, bs = "ps"),
         ~ s(sqrt_cape, bs = "ps") + s(hcc, bs = "ps") + s(str, bs = "ps") +
           s(tciw, bs = "ps") + s(q_prof_PC1, bs = "ps")
)

## Estimate models.
set.seed(123)
flash_model_ztnbinom <- bamlss(f, data = FlashAustriaTrain, ## Standard interface.
  family = "ztnbinom", binning = TRUE,                      ## General arguments.
  optimizer = opt_boost, maxit = 1000,                      ## Boosting arguments.
  thin = 10, burnin = 2000, n.iter = 12000,                 ## Sampler arguments.
  light = TRUE)

set.seed(123)
flash_model_ztSICHEL <- bamlss(f, data = FlashAustriaTrain,
  family = ztSICHEL, binning = TRUE,
  optimizer = opt_boost, maxit = 1000,
  thin = 50, burnin = 2000, n.iter = 3000, cores = 50,
  light = TRUE)

set.seed(123)
flash_model_ztBNB <- bamlss(f, data = FlashAustriaTrain,
  family = ztBNB, binning = TRUE,
  optimizer = opt_boost, maxit = 1000,
  thin = 50, burnin = 2000, n.iter = 3000, cores = 50,
  light = TRUE)

## Reduced formula for ztnbinom model based on boosting step.
boost_summary(flash_model_ztnbinom)

fr <- list(
  counts ~ s(sqrt_cape, bs = "ps") + s(hcc, bs = "ps") + s(str, bs = "ps") +
           s(tciw, bs = "ps") + s(q_prof_PC1, bs = "ps"),
         ~ s(sqrt_cape, bs = "ps")
)

set.seed(123)
flash_model_ztnbinom_r <- bamlss(fr, data = FlashAustriaTrain,
  family = "ztnbinom", binning = TRUE,
  optimizer = opt_boost, maxit = 1000,
  thin = 10, burnin = 2000, n.iter = 12000,
  light = TRUE)

save(
  flash_model_ztnbinom_r,
  flash_model_ztnbinom,
  flash_model_ztSICHEL,
  flash_model_ztBNB,
  file = "FlashAustriaModel.rda"
)

