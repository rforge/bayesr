## source all functions
sbayesr()

## generate data
n <- 500
dat <- data.frame("x1" = runif(n, -3, 3))
dat$fac <- factor(rep(1:50, length.out = n))
dat$x2 <- with(dat, runif(length(unique(fac)), -3, 3)[fac])
dat$re <- with(dat, cos(x2) + rnorm(length(unique(fac)), sd = 0.1)[fac])

fun <- function(x, theta = c(2, -20, -0.1)) {
  theta[1] * exp(theta[2] * exp(theta[3] * x))
}

## dat$y <- with(dat, 21 + cf + sin(x2) + rnorm(n, sd = 0.1))
dat$y <- hs(with(dat, 1.2 + sin(x1) + re + rnorm(n, sd = 0.1)), 0.0001, 2)

## fit model
b <- bayesx2(y ~ sx(x1) + sx(fac, bs = "re") & fac ~ sx(x2) | y ~ sx(fac, bs = "re") & fac ~ sx(x2), family = normal.BayesX, data = dat, dir = "~/tmp")


b <- bayesr(y ~ s(x1) + s(fac, bs = "re") | fac ~ s(x2), data = dat)

dat <- dat[order(dat$x1), ]

X <- mgcv:::smooth2random(smoothCon(s(x1), dat, NULL, absorb.cons = TRUE)[[1]], names(dat), type = 2)
Xf <- X$Xf
Xr <- X$rand$Xr

dat2 <- data.frame("fixed" = Xf, "random" = Xr)
dat2$y <- dat$y

b <- bayesx(y ~ fixed + sx(random.1, bs = "ridge") + sx(random.2, bs = "ridge") +
  sx(random.3, bs = "ridge") + sx(random.4, bs = "ridge") + sx(random.5, bs = "ridge") +
  sx(random.6, bs = "ridge") + sx(random.7, bs = "ridge") + sx(random.8, bs = "ridge"),
  data = dat2)



p <- predict(b, term = "s(x1)")
plot(dat$x1, p)
points(dat$x1, dat$cf * fun(dat$x1), col = 2, pch = 3)


## generate data
hs <- function(x, min = 0.1, max = 0.6) {
  x <- if(length(unique(x)) > 1) {
    (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE)) * (max - min) + min
  } else x
  x
}

n <- 500
dat <- data.frame("x1" = runif(n, -3, 3), x2 = runif(n, -3, 3))
dat$y <- hs(with(dat, 1.2 + sin(x1) + cos(x2) + rnorm(n, sd = (cos(dat$x1) + 2) / 4)), 0.001, 0.999)

a <- bayesr(y ~ s(x1) + s(x2) | s(x1), data = dat, family = gaussian.JAGS)

b <- bayesx2(y ~ sx(x1) + sx(x2) | sx(x1), data = dat, family = normal)

