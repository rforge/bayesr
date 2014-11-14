###########
## Bug 1 ##
###########
## Main text when plotting 2nd stage hierarchical models
n <- 500
dat <- data.frame("x1" = runif(n, -3, 3))
dat$fac <- factor(rep(1:50, length.out = n))
dat$fac2 <- factor(rep(1:5, length.out = n))
dat$x2 <- with(dat, runif(length(unique(fac)), -3, 3)[fac])
dat$re <- with(dat, cos(x2) + rnorm(length(unique(fac)), sd = 0.1)[fac])
dat <- dat[order(as.integer(dat$fac2)), ]
dat$y <- scale2(with(dat, 1.2 + sin(x1) + re + rnorm(n, sd = 0.1)), 0.0001, 2)

f <- list(
  y ~ -1 + sx(x1) + sx(fac, bs = "re") + sx(fac2, bs = "re"),
  fac ~ 1 + sx(x2),
  fac2 ~ -1 + sx(x2),
  sigma2 ~ sx(x1)
)

b <- bamlss(f, family = gaussian2, data = dat, engine = "BayesX")
plot(b)
