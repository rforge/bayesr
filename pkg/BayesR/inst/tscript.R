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
dat$y <- with(dat, 1.2 + sin(x1) + re + rnorm(n, sd = 0.1))

## fit model
b <- bayesx2(y ~ s(x1) + s(fac, bs = "re") | fac ~ s(x2), data = dat)
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
