## generate data
set.seed(123)
n <- 500
dat <- data.frame("x1" = runif(n, -3, 3))
dat$fac <- factor(rep(1:10, length.out = n))
rc <- rnorm(nlevels(dat$fac), sd = 0.6) + 1
dat$y <- with(dat, 1.2 + sin(x1) * rc[fac] + rnorm(n, sd = 0.1))


## fit model
f <- list(
  y ~ sx(x1, bs = "rps", by = fac, sum2 = 10),
  fac ~ -1
)

b <- bayesr(f, family = gaussian2, data = dat, engine = "BayesX")

plot(b, residuals = TRUE)



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

a <- bayesr(y ~ s(x1) + s(x2), ~ s(x1), data = dat, method = "MCMC")

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
f <- simfun(type = "pick")

n <- 100
dat <- data.frame("x1" = sort(runif(n, 0, 1)))
dat$y <- with(dat, 1.2 + f(x1) + rnorm(n, sd = 0.1))

b <- bayesr(y ~ -1 + rs(x1, bs = "ps", fx = TRUE, xt = list(center = FALSE)), data = dat, engine = "JAGS")



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
dat$y <- with(dat, 1.2 + sin(x1) * cos(x2) + rnorm(n, sd = scale2(cos(x1), 0.1, 0.8)))

b <- bayesr(y ~ s(x1, x2), ~ s(x1), data = dat, method = "MCMC", svalues = FALSE, n.iter = 1200, burnin = 200, thin = 1, family = gaussian2, propose = "slice")

plot(b)
plot(b, which = "samples")


## random test
n <- 500
dat <- data.frame("x" = sort(runif(n, -3, 3)))
dat$y <- with(dat, 1.2 + sin(x) + rnorm(n, sd = scale3(cos(x), 0.1, 0.8)))

X <- mgcv:::smooth2random(smoothCon(s(x), dat, NULL, absorb.cons = TRUE)[[1]], names(dat), type = 2)
Xf <- X$Xf
Xr <- X$rand$Xr



## Testing numerical derivatives.
library("BayesR")
library("gamlss")

n <- 500
x <- seq(-3, 3, length = n)
mu <- sin(x)
sigma <- scale2(x^2, 0.2, 0.5)
y <- scale2(rnorm(n, mu, sigma), 0.01, 0.99)

family <- tF(NO)

eta <- list(mu = mu, sigma = sigma)

par(mfrow = c(2, 2))

plot(num_deriv(y, eta, family, id = "mu") ~ family$score$mu(y, eta))
plot(num_deriv(y, eta, family, id = "sigma") ~ family$score$sigma(y, eta))

plot(-1 * num_deriv(y, eta, family, id = "mu", d = 2) ~ family$weights$mu(y, eta))
plot(-1 * num_deriv(y, eta, family, id = "sigma", d = 2) ~ family$weights$sigma(y, eta))


