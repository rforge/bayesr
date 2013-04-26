## source all functions
sbayesr()

## generate data
n <- 500
dat <- data.frame("x1" = runif(n, -3, 3), "x2" = runif(n, -3, 1))
dat$fac <- factor(rep(1:5, length.out = n))
dat$cf <- with(dat, rnorm(5, sd = 1.2)[fac])

fun <- function(x, theta = c(2, -20, -0.1)) {
  theta[1] * exp(theta[2] * exp(theta[3] * x))
}

## dat$y <- with(dat, 21 + cf + sin(x2) + rnorm(n, sd = 0.1))
dat$y <- with(dat, 1.2 + sin(x1) + rnorm(n, sd = 0.1))

## fit model
b <- bayesr(y ~ s(x1, k = 10, bs = "ps"), data = dat)


p <- predict(b, term = "s(x1)")
plot(dat$x1, p)
points(dat$x1, dat$cf * fun(dat$x1), col = 2, pch = 3)
