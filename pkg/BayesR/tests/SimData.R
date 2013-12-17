library("BayesR")

## Create a grid map.
gm <- pixelmap(n = 15)
gmap <- gm$map

## Create the data set.
set.seed(111)
dat <- gm$data[rep(seq_along(gm$data[, 1]), 5), ]
n <- nrow(dat)
dat$x1 <- runif(n, -3, 3)
dat$x2 <- runif(n, -3, 3)
dat$sp <- with(dat, sin(x) * cos(y))
dat$sp <- dat$sp - mean(dat$sp)
dat$re <- with(dat, rnorm(length(unique(id)), sd = 0.6)[id])

## The linear predictor.
dat$eta <- with(dat, 1.5 + sin(x1) + cos(x2) + sp + re)

## Gaussian response.
dat$yg <- with(dat, eta + rnorm(n, sd = 0.6))

## Model with JAGS.
b0 <- bayesr(yg ~ s(x1) + s(x2) + s(id, bs = "mrf", xt = list(penalty = gm$nmat)) +
  s(id, bs = "re"), data = dat, family = gaussian, engine = "JAGS")

## Model with BayesX.
b1 <- bayesr(yg ~ sx(x1) + sx(x2) + sx(id, bs = "mrf", map = gm$map) +
  sx(id, bs = "re"), id ~ -1, data = dat, family = gaussian, engine = "BayesX")
