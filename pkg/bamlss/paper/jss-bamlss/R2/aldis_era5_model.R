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
set.seed(456)
b_ztnbinom <- bamlss(f, data = FlashAustriaTrain,       ## Standard interface.
  family = "ztnbinom", binning = TRUE,                  ## General arguments.
  optimizer = opt_boost, maxit = 1000,                  ## Boosting arguments.
  thin = 50, burnin = 500, n.iter = 1500, cores = 50,   ## Sampler arguments.
  light = TRUE)

set.seed(456)
b_ztSICHEL <- bamlss(f, data = FlashAustriaTrain,
  family = ztSICHEL, binning = TRUE,
  optimizer = opt_boost, maxit = 1000,
  thin = 50, burnin = 500, n.iter = 1500, cores = 50,
  light = TRUE)

set.seed(456)
b_ztBNB <- bamlss(f, data = FlashAustriaTrain,
  family = ztBNB, binning = TRUE,
  optimizer = opt_boost, maxit = 1000,
  thin = 50, burnin = 500, n.iter = 1500, cores = 50,
  light = TRUE)

save(b_ztnbinom, b_ztSICHEL, b_ztBNB, file = "flashmodels.rda")

### --- predictions ---
fit1 <- predict(b_ztnbinom, newdata = FlashAustriaEval, type = "parameter")
fit2 <- predict(b_ztSI, newdata = FlashAustriaEval, type = "parameter")
fit3 <- predict(b_ztBNB, newdata = FlashAustriaEval, type = "parameter")

family(b_ztnbinom)$loglik(FlashAustriaEval$counts, fit1)
family(b_ztSI)$loglik(FlashAustriaEval$counts, fit2)
family(b_ztBNB)$loglik(FlashAustriaEval$counts, fit3)

e1t <- residuals(b_ztnbinom, newdata = FlashAustriaTrain)
e2t <- residuals(b_ztSI, newdata = FlashAustriaTrain)
e3t <- residuals(b_ztBNB, newdata = FlashAustriaTrain)

e1e <- residuals(b_ztnbinom, newdata = FlashAustriaEval)
e2e <- residuals(b_ztSI, newdata = FlashAustriaEval)
e3e <- residuals(b_ztBNB, newdata = FlashAustriaEval)

par(mfrow = c(1, 2))
plot(c("SItr" = e1t, "ztnbinom" = e2t, "BNB" = e3t), which = "wp", main = "Training", pos = "top")
plot(c("SItr" = e1e, "ztnbinom" = e2e, "BNB" = e3e), which = "wp", main = "Testing", pos = "top")

