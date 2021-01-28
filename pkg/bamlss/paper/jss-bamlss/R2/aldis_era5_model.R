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

## Model formula.
f <- list(
  counts ~ s(d2m, bs = "ps") + s(q_prof_PC1, bs = "ps") +
           s(cswc_prof_PC4, bs = "ps") + s(t_prof_PC1, bs = "ps") +
           s(v_prof_PC2, bs = "ps") + s(sqrt_cape, bs = "ps"), 
         ~ s(sqrt_lsp, bs = "ps")
         ~ 1
)

## Estimate models.
set.seed(123)
flash_model_ztnbinom <- bamlss(f, data = FlashAustriaTrain, ## Standard interface.
  family = "ztnbinom", binning = TRUE,                      ## General arguments.
  optimizer = opt_boost, maxit = 1000,                      ## Boosting arguments.
  thin = 3, burnin = 1000, n.iter = 3000,                 ## Sampler arguments.
  light = TRUE)

set.seed(123)
flash_model_ztSICHEL <- bamlss(f, data = FlashAustriaTrain,
  family = ztSICHEL, binning = TRUE,
  optimizer = opt_boost, maxit = 1000,
  thin = 3, burnin = 1000, n.iter = 3000,
  light = TRUE)

save(
  flash_model_ztnbinom,
  flash_model_ztSICHEL,
  file = "FlashAustriaModel.rda"
)

