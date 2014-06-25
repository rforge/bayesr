## generate data
set.seed(123)
n <- 500
dat <- data.frame("x1" = runif(n, -3, 3))
dat$fac <- factor(rep(1:10, length.out = n))
rc <- rnorm(nlevels(dat$fac), sd = 0.3) + 1
dat$y <- with(dat, 1.2 + x1 * rc[fac] + rnorm(n, sd = 0.1))


## fit model
f <- list(
  y ~ -1 + sx(fac, by = x1, bs = "re"),
  fac ~ 1
)

b <- bayesr(f, family = gaussian2, data = dat, engine = "BayesX")

plot(b)



n <- 50
nd <- data.frame(x1 = seq(-10, 10, length = n), x2 = seq(-10, 10, length = n))
nd$p <- predict(b, nd, model = c("mu", "h1"), term = "sx(x1)")

nd$p <- predict(b2, nd, model = "mu", term = "sx(x1)")


b <- bayesr(y ~ s(x1) + s(fac, bs = "re") | fac ~ s(x2), data = dat)

dat <- dat[order(dat$x1), ]

X <- mgcv:::smooth2random(smoothCon(s(x1), dat, NULL, absorb.cons = TRUE)[[1]], names(dat), type = 2)
Xf <- X$Xf
Xr <- X$rand$Xr

dat2 <- data.frame("fixed" = Xf, "random" = Xr)
dat2$y <- dat$y

b <- bayesx2(y ~ fixed + sx(random.1, bs = "ridge") + sx(random.2, bs = "ridge") +
  sx(random.3, bs = "ridge") + sx(random.4, bs = "ridge") + sx(random.5, bs = "ridge") +
  sx(random.6, bs = "ridge") + sx(random.7, bs = "ridge") + sx(random.8, bs = "ridge"),
  data = dat2)


p <- predict(b, term = "s(x1)")
plot(dat$x1, p)
points(dat$x1, dat$cf * fun(dat$x1), col = 2, pch = 3)


## generate data
library("gamlss")

set.seed(111)
n <- 200
dat <- data.frame("x1" = sort(runif(n, -3, 3)), x2 = runif(n, -3, 3))
dat$y <- scale2(with(dat, 1.2 + sin(x1) + cos(x2) + rnorm(n, sd = (cos(dat$x1) + 2) / 4)), 0.001, 0.999)

a <- bayesr(y ~ s(x1, x2, bs = "kr", k = 50), ~ s(x1), data = dat, update = "iwls", propose = "iwls", method = c("backfitting", "MCMC"))

plot(a)


b <- bayesr(y ~ s(x1) + s(x2), ~ s(x1), data = dat, family = tF(BE), method = "MCMC")

plot(c(a, b))


a <- bayesr(y1 | y2 ~ s(x1) + s(x2), data = dat, family = gaussian.BayesR)

b2 <- bayesx2(y ~ sx(x1) + sx(x2), ~ sx(x1), data = dat)


data("GAMart")
b <- bayesr(cat ~ s(x1) + s(x2) + s(x3) + s(long, lat), family = multinomial.BayesR, data = GAMart)



data("marital.nz", package = "VGAM")
b <- bayesr(mstatus ~ s(age), family = multinomial.BayesR, data = marital.nz)



## pick function
f <- simfun(type = "complicated")

n <- 200
dat <- data.frame("x1" = sort(runif(n, 0, 1)))
dat$y <- with(dat, 1.2 + f(x1) + rnorm(n, sd = 0.2))

b <- bayesr(y ~ rs(x1,fx=TRUE), data = dat, method = "backfitting")

g <- coef(b)
g <- g[grep("s(x1)", names(g), fixed = TRUE)]

X <- smooth.construct(rs(x1), dat, NULL)

dat$f <- X$get.mu(X$X, g)

dat$y2 <- with(dat, 1.2 + f + rnorm(n, sd = 0.5))

b2 <- bayesr(y2 ~ rs(x1,fx=TRUE), data = dat, method = "backfitting")

g2 <- coef(b2)
g2 <- g2[grep("s(x1)", names(g2), fixed = TRUE)]

plot(g, g2)
abline(a = 0, b = 1)


f <- simfun(type = "2d")
n <- 200
dat <- data.frame("x1" = sort(runif(n, 0, 1)), "x2" = runif(n, 0, 1))
dat$y <- with(dat, 1.2 + f(x1, x2) + rnorm(n, sd = 0.2))

b0 <- bayesr(y ~ s(x1, k = 6,fx=T):s(x2, k = 6,fx=T), data = dat, method = c("backfitting2", "MCMC"))
b1 <- bayesr(y ~ s(x1, x2, k = 20), data = dat)

plot(c(b0, b1), type = "mba")



library("MASS")
data("mcycle")
mcycle$accel2 <- scale2(mcycle$accel, -1, 1)

b <- bayesr(accel2 ~ rs(times, k = 20), ~ s(times, bs = "ps"), data = mcycle, method = "backfitting2")


## by variable test
n <- 500
dat <- data.frame(x1 = sort(runif(n, -3, 3)), x2 = sort(runif(n)))
dat$y <- with(dat, 1.2 + sin(x1) * x2 + rnorm(n, sd = 0.2))

b0 <- gam(y ~ x2 + s(x1, by = x2), data = dat)
b1 <- bayesr(y ~ x2 + s(x1, by = x2), data = dat, method = "MCMC")
b2 <- bayesr(y ~ x2 + sx(x1, by = x2), data = dat, engine = "BayesX")


## IWLS test
n <- 100
dat <- data.frame("x1" = runif(n, -3, 3), "x2" = runif(n, -3, 3))
dat$y <- with(dat, 1.2 + sin(x1) + rnorm(n, sd = scale2(cos(x1), 0.1, 0.8)))

b0 <- bayesr(y ~ s(x1), ~ s(x1), data = dat, method = c("backfitting", "MCMC"),
  n.iter = 12000, burnin = 2000, thin = 10, family = gaussian2,
  update = "iwls", propose = "wslice")

b1 <- bayesr(y ~ s(x1), ~ s(x1), data = dat, method = c("backfitting", "MCMC"),
  n.iter = 12000, burnin = 2000, thin = 10, family = gaussian2,
  update = "optim", propose = "oslice")

b2 <- bayesr(y ~ s(x1), ~ s(x1), data = dat, method = c("backfitting", "MCMC"),
  n.iter = 12000, burnin = 2000, thin = 10, family = gaussian2,
  update = "optim", propose = "slice")

b3 <- bayesr(y ~ s(x1), ~ s(x1), data = dat, method = c("backfitting", "MCMC"),
  family = gaussian2, update = "iwls", propose = "iwls")

plot(b)
plot(b, which = "samples")


## random test
n <- 100
dat <- data.frame("x" = sort(runif(n, -3, 3)))
dat$y <- with(dat, 1.2 + sin(x) + rnorm(n, sd = scale3(cos(x), 0.1, 0.8)))

X <- mgcv:::smooth2random(smoothCon(s(x), dat, NULL, absorb.cons = TRUE)[[1]], names(dat), type = 2)
Xf <- X$Xf
Xr <- X$rand$Xr
U <- X$trans.U
D <- X$trans.D



## Testing numerical derivatives.
library("BayesR")
library("gamlss")

n <- 50
x <- seq(-3, 3, length = n)
mu <- scale2(sin(x), 0.5, 1.5)
sigma <- scale2(x^2, 0.2, 0.5)
y <- scale2(rnorm(n, mu, sigma), 0.01, 0.99)

family <- tF(NO2)

eta <- list(mu = mu, sigma = sigma)

par(mfrow = c(2, 2))

plot(num_deriv(y, eta, family, id = "mu") ~ family$iwls$score$mu(y, eta))
points(num_deriv2(y, eta, family, id = "mu") ~ family$iwls$score$mu(y, eta), col = 2)
plot(num_deriv(y, eta, family, id = "sigma") ~ family$iwls$score$sigma(y, eta))
points(num_deriv2(y, eta, family, id = "sigma") ~ family$iwls$score$sigma(y, eta), col = 2)

plot(num_deriv(y, eta, family, id = "mu", d = 2) ~ family$iwls$weights$mu(y, eta))
points(num_deriv2(y, eta, family, id = "mu", d = 2) ~ family$iwls$weights$mu(y, eta), col = 2)
plot(num_deriv(y, eta, family, id = "sigma", d = 2) ~ family$iwls$weights$sigma(y, eta))
points(num_deriv2(y, eta, family, id = "sigma", d = 2) ~ family$iwls$weights$sigma(y, eta), col = 2)



## More tests.
d <- dgp_gaussian2(
  300,
  mu = list(const = -0.5, type = list(c("double", "spatial", "const"))),
  sigma2 = list(const = 0.01, type = list("quadratic", "const"))
)

f <- list(
  y ~ s(mu.x11) + s(mu.long1, mu.lat1, k = 50, bs = "kr"),
  sigma2 ~ s(sigma2.x11)
)

b0 <- bayesr(f, data = d, update = "iwls", propose = "iwls", method = c("backfitting", "MCMC"))
b1 <- bayesr(f, data = d, update = "iwls", propose = "wslice", method = c("backfitting", "MCMC"))

f <- list(
  y ~ sx(mu.x11) + sx(mu.long1, mu.lat1, knots = 50, bs = "kr", update = "orthogonal"),
  sigma2 ~ sx(sigma2.x11)
)

b2 <- bayesr(f, data = d, engine = "BayesX", family = gaussian2)



d <- dgp_beta(
  200,
  mu = list(const = -0.5, type = list(c("pick", "const"))),
  sigma2 = list(const = 0.01, type = list("quadratic", "const"))
)

no <- gaussian2.BayesR()
no$iwls <- NULL

f <- list(
  y ~ s(mu.x11),
  sigma ~ s(sigma2.x11)
)

b <- bayesr(f, data = d, family = tF(BE), update = "iwls")


