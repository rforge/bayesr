## R2BayesX testing
library("R2BayesX")


## (1) model including factors
set.seed(111)
n <- 200
     
## regressors
dat <- data.frame(x = runif(n, -3, 3), fac = factor(rep(1:10, n/10), labels = letters[1:10]))
     
## response
dat$y <- with(dat, 1.5 + sin(x) + c(2.67, 5, 6, 3, 4, 2, 6, 7, 9, 7.5)[fac] + rnorm(n, sd = 0.6))
     
## estimate models
b1 <- bayesx(y ~ sx(x) + fac, data = dat, method = "MCMC")
b2 <- bayesx(y ~ sx(x) + fac, data = dat, method = "REML")
b3 <- bayesx(y ~ sx(x) + fac, data = dat, method = "STEP")

## summaries
summary(b1)
summary(b2)
summary(b3)


## (2) models with transformations
b1 <- bayesx(log(y) ~ sx(scale(x)) + fac, data = dat, method = "MCMC")
b2 <- bayesx(log(y) ~ sx(scale(x)) + fac, data = dat, method = "REML")
b3 <- bayesx(log(y) ~ sx(scale(x)) + fac, data = dat, method = "STEP")

## summaries
summary(b1)
summary(b2)
summary(b3)

