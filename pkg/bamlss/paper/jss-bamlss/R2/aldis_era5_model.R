## Set the seed for reproducibility. 
set.seed(123)

## Required packages.
library("bamlss")
library("gamlss.dist")
library("gamlss.tr")

## Load the data.
data("FlashAustria", package = "FlashAustria")

## Generate zero truncated Sichel distribution.
ztSI <- trun(0, family = "SI", local = FALSE)

## Generate zero truncated beta negative binomial.
ztBNB <- trun(0, family = "BNB", local = FALSE)

## Model formula, up to four parameters.
f <- list(
  counts ~ s(sqrt_cape) + s(hcc) + s(str) + s(tciw) + s(q_prof_PC1),
         ~ s(sqrt_cape) + s(hcc) + s(str) + s(tciw) + s(q_prof_PC1),
         ~ s(sqrt_cape) + s(hcc) + s(str) + s(tciw) + s(q_prof_PC1),
         ~ s(sqrt_cape) + s(hcc) + s(str) + s(tciw) + s(q_prof_PC1)
)

## Estimate models.
b_ztnbinom <- bamlss(f, data = FlashAustriaTrain, ## Standard interface.
  family = "ztnbinom", binning = TRUE,            ## General arguments.
  optimizer = opt_boost, maxit = 1000,            ## Boosting arguments.
  thin = 5, burnin = 1000, n.iter = 6000,         ## Sampler arguments.
  light = TRUE)

b_ztSI <- bamlss(f, data = FlashAustriaTrain,
  family = ztSI, binning = TRUE,
  optimizer = opt_boost, maxit = 1000,
  thin = 5, burnin = 1000, n.iter = 6000,
  light = TRUE)

b_ztBNB <- bamlss(f, data = FlashAustriaTrain,
  family = ztBNB, binning = TRUE,
  optimizer = opt_boost, maxit = 1000,
  thin = 5, burnin = 1000, n.iter = 6000,
  light = TRUE)

save(b_ztnbinom, b_ztSI, b_ztBNB, file = "flashmodels.rda")

### --- predictions ---
#fit1 <- predict(b1, newdata = d_eval, type = "parameter")
#fit2 <- predict(b2, newdata = d_eval, type = "parameter")

#family(b1)$loglik(d_eval$counts, fit1)
#family(b2)$loglik(d_eval$counts, fit2)

#e1t <- residuals(b1, newdata = FlashAustriaTrain)
#e2t <- residuals(b2, newdata = FlashAustriaTrain)

#e1e <- residuals(b1, newdata = FlashAustriaEval)
#e2e <- residuals(b2, newdata = FlashAustriaEval)

#par(mfrow = c(1, 2))
#plot(c("SItr" = e1t, "ztnbinom" = e2t), which = "wp", main = "Training", pos = "top")
#plot(c("SItr" = e1e, "ztnbinom" = e2e), which = "wp", main = "Testing", pos = "top")

