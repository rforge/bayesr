## generate data
n <- 500
dat <- data.frame("x1" = runif(n, -3, 3))
dat$fac <- factor(rep(1:50, length.out = n))
dat$fac2 <- factor(rep(1:5, length.out = n))
dat$x2 <- with(dat, runif(length(unique(fac)), -3, 3)[fac])
dat$re <- with(dat, cos(x2) + rnorm(length(unique(fac)), sd = 0.1)[fac])

fun <- function(x, theta = c(2, -20, -0.1)) {
  theta[1] * exp(theta[2] * exp(theta[3] * x))
}

## dat$y <- with(dat, 21 + cf + sin(x2) + rnorm(n, sd = 0.1))
dat$y <- scale2(with(dat, 1.2 + sin(x1) + re + rnorm(n, sd = 0.1)), 0.0001, 2)

## fit model
b <- bayesx2(y ~ x1 + fac, family = gaussian.BayesR, data = dat)

b <- bayesx2(y ~ x1 + sx(fac, bs = "re") + fac2, fac ~ - 1 + sx(x2), ~ sx(x1) + sx(fac, bs = "re"), family = gaussian.BayesR, data = dat)

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
n <- 500
dat <- data.frame("x1" = sort(runif(n, -3, 3)), x2 = runif(n, -3, 3))
dat$y <- scale2(with(dat, 1.2 + sin(x1) + cos(x2) + rnorm(n, sd = (cos(dat$x1) + 2) / 4)), 0.001, 0.999)

a <- bayesr(y ~ s(x1) + s(x2), ~ s(x1), data = dat, engine = "IWLS", family = beta)


a <- bayesr(y1 | y2 ~ s(x1) + s(x2), data = dat, family = gaussian.BayesR)

b2 <- bayesx2(y ~ sx(x1) + sx(x2), ~ sx(x1), data = dat)


data("GAMart")
b <- bayesr(cat ~ s(x1) + s(x2) + s(x3) + s(long, lat), family = multinomial.BayesR, data = GAMart)



data("marital.nz", package = "VGAM")
b <- bayesr(mstatus ~ s(age), family = multinomial.BayesR, data = marital.nz)



## pick function
funpick <- function(x, min = 0, max = 1) { 
  x <- (x - min) / (max - min)
  y <- sin(2 * (4 * x - 2)) + 2 * exp(-16^2 * (x - 0.5)^2)
  return(y - mean(y))
}

n <- 500
dat <- data.frame("x1" = sort(runif(n, 0, 1)))
dat$y <- with(dat, 1.2 + funpick(x1) + rnorm(n, sd = 0.1))

b <- bayesr(y ~ -1 + rs(x1, bs = "ps", k = 5, fx = TRUE, xt = list(center = FALSE)), data = dat,
  sampler = list(n.iter = 22000, burnin = 12000, thin = 2))



## by variable test
n <- 500
dat <- data.frame(x1 = sort(runif(n, -3, 3)), x2 = sort(runif(n)))
dat$y <- with(dat, 1.2 + sin(x1) * x2 + rnorm(n, sd = 0.2))

b0 <- gam(y ~ s(x1, by = x2), data = dat)
b1 <- bayesr(y ~ s(x1, by = x2), data = dat)
b2 <- bayesx2(y ~ sx(x1, by = x2), data = dat)


## IWLS test
n <- 500
dat <- data.frame("x" = sort(runif(n, -3, 3)))
dat$y <- with(dat, 1.2 + sin(x) + rnorm(n, sd = scale3(cos(x), 0.1, 0.8)))

b <- bayesr(y ~ s(x), ~ s(x), family = gaussian.BayesR, data = dat)

plot(b)
plot(b, which = "samples")

