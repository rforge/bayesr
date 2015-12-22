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

b <- bamlss(f, family = gaussian2, data = dat, engine = "BayesX")

plot(b)



n <- 50
nd <- data.frame(x1 = seq(-10, 10, length = n), x2 = seq(-10, 10, length = n))
nd$p <- predict(b, nd, model = c("mu", "h1"), term = "sx(x1)")

nd$p <- predict(b2, nd, model = "mu", term = "sx(x1)")


b <- bamlss(y ~ s(x1) + s(fac, bs = "re") | fac ~ s(x2), data = dat)

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

a <- bamlss(y ~ s(x1, x2, bs = "kr", k = 50), ~ s(x1), data = dat, update = "iwls", propose = "iwls", method = c("backfitting", "MCMC"))

plot(a)


b <- bamlss(y ~ s(x1) + s(x2), ~ s(x1), data = dat, family = tF(BE), method = "MCMC")

plot(c(a, b))


a <- bamlss(y1 | y2 ~ s(x1) + s(x2), data = dat, family = gaussian.bamlss)

b2 <- bayesx2(y ~ sx(x1) + sx(x2), ~ sx(x1), data = dat)


data("GAMart")
b <- bamlss(cat ~ s(x1) + s(x2) + s(x3) + s(long, lat), family = multinomial.bamlss, data = GAMart)



data("marital.nz", package = "VGAM")
b <- bamlss(mstatus ~ s(age), family = multinomial.bamlss, data = marital.nz)



## pick function
f <- simfun(type = "pick")

set.seed(111)
n <- 200
dat <- data.frame("x1" = sort(runif(n, 0, 1)))
dat$y <- with(dat, 10 + f(x1) + rnorm(n, sd = 0.2))

k0 <- 8
X <- smoothCon(s(x1, k = k0, bs = "tp"), dat, NULL, absorb.cons = TRUE)[[1]]$X
k <- ncol(X)

ff <- function(beta) {
  f <- (X %*% beta[2:(k + 1)]) / exp(1 + X %*% beta[(k + 2):(2 * k + 1)])
  f <- f - mean(f)
  beta[1] + f
}

objfun <- function(beta) {
  f <- ff(beta)
  return(sum((dat$y - f)^2))
}

start <- c(0, rep(0, ncol(X) * 2))

opt <- optim(start, objfun, method = "BFGS")

plot(dat)
lines(ff(opt$par) ~ dat$x1)
b <- gam(y ~ s(x1, k = k0), data = dat)
lines(fitted(b) ~ dat$x1, col = "red")


X2 <- cbind(X, dat$y * X)
b <- lm(y ~ X2, data = dat)
lines(ff(coef(b)) ~ dat$x1, col = "blue")


stepfun <- function(maxit = 100, eps = .Machine$double.eps)
{
  y <- dat$y
  k <- ncol(X)

  beta <- rep(eps, 2 * k + 1)

  fu <- drop(beta[1] + X %*% beta[2:(k + 1)])
  fl <- drop(1 + X %*% beta[(k + 2):(2 * k + 1)])
  fit <- fu / fl

  iter <- 0; err <- eps + 1
  while(iter < maxit & err > eps) {
    fit0 <- fit

    bu <- coef(lm.fit(cbind(1, X) / fl, y))
    fu <- drop(bu[1] + X %*% bu[2:(k + 1)])
    fit <- fu / fl

    bl <- coef(lm.fit(X, y * fl))
    fl <- drop(1 + X %*% bl)
    fit <- fu / fl

    if(iter > 0)
      err <- mean((fit - fit0) / fit0)

    iter <- iter + 1
  }

  b <- c(bu, bl)
  plot(dat, ylim = range(c(fit, dat$y)))
  lines(fit ~ dat$x1, col = "green")
}

stepfun()








a <- nls(y ~ ns(x1,df=5) / ns(x1,df=5), data = dat)


b <- bamlss(y ~ s(x1), data = dat, method = "MP")

b <- bamlss(y ~ s(x1), data = dat, method = "MCMC", sample = "slice", svalues = FALSE, n.iter = 1200, burnin = 200, thin = 1)

g <- coef(b)
g <- g[grep("s(x1)", names(g), fixed = TRUE)]

X <- smooth.construct(rs(x1), dat, NULL)

dat$f <- X$get.mu(X$X, g)

dat$y2 <- with(dat, 1.2 + f + rnorm(n, sd = 0.2))

b2 <- bamlss(y2 ~ rs(x1), data = dat, method = "backfitting")

g2 <- coef(b2)
g2 <- g2[grep("s(x1)", names(g2), fixed = TRUE)]

plot(g, g2)
abline(a = 0, b = 1)


f <- simfun(type = "2d")
n <- 500
dat <- data.frame("x1" = sort(runif(n, 0, 1)), "x2" = runif(n, 0, 1))
dat$y <- with(dat, 1.2 + f(x1, x2) + rnorm(n, sd = 0.2))

b0 <- bamlss(y ~ rs(s(x1), s(x2), link = "inverse"), data = dat, method = "backfitting")
b1 <- bamlss(y ~ s(x1, x2, k = 20), data = dat, method = "backfitting")

plot(c(b0, b1), type = "mba")



library("MASS")
data("mcycle")
mcycle$accel2 <- scale(mcycle$accel)

b <- bamlss(accel2 ~ s(times, k = 20), data = mcycle, method = "backfitting")


## by variable test
n <- 500
dat <- data.frame(x1 = runif(n, -3, 3), x2 = runif(n),
  fac = as.factor(sample(1:2, n, replace = TRUE)))
dat$y <- with(dat, 1.2 + sin(x1) * x2 + rnorm(n, sd = 0.2))

X <- smooth.construct(s(x1, by = x2), dat, NULL)$X
b <- solve(crossprod(X)) %*% t(X) %*% dat$y
f <- X %*% b
plot2d(f ~ dat$x1)


b0 <- bamlss(y ~ s(x1, by = x2), data = dat)
b1 <- bamlss(y ~ x2 + s(x1, by = x2), data = dat, method = "MCMC")
b2 <- bamlss(y ~ x2 + sx(x1, by = x2), data = dat, engine = "BayesX", family = gaussian)


## IWLS test
n <- 500
dat <- data.frame("x1" = runif(n, -3, 3), "x2" = runif(n, -3, 3))
dat$y <- with(dat, 1.2 + sin(x1) + rnorm(n, sd = scale2(cos(x1), 0.1, 0.8)))

b0 <- bamlss(y ~ s(x1), ~ s(x1), data = dat, family = gaussian2, method = "MP2")

b0 <- bamlss(y ~ s(x1), ~ s(x1), data = dat, method = c("backfitting", "MCMC"),
  n.iter = 12000, burnin = 2000, thin = 10, family = gaussian2,
  update = "iwls", propose = "wslice")

b1 <- bamlss(y ~ s(x1), ~ s(x1), data = dat, method = c("backfitting", "MCMC"),
  n.iter = 12000, burnin = 2000, thin = 10, family = gaussian2,
  update = "optim", propose = "oslice")

b2 <- bamlss(y ~ s(x1), ~ s(x1), data = dat, method = c("backfitting", "MCMC"),
  n.iter = 12000, burnin = 2000, thin = 10, family = gaussian2,
  update = "optim", propose = "slice")

b3 <- bamlss(y ~ s(x1), ~ s(x1), data = dat, method = c("backfitting", "MCMC"),
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
library("bamlss")
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

b0 <- bamlss(f, data = d, update = "iwls", propose = "iwls", method = c("backfitting", "MCMC"))
b1 <- bamlss(f, data = d, update = "iwls", propose = "wslice", method = c("backfitting", "MCMC"))

f <- list(
  y ~ sx(mu.x11) + sx(mu.long1, mu.lat1, knots = 50, bs = "kr", update = "orthogonal"),
  sigma2 ~ sx(sigma2.x11)
)

b2 <- bamlss(f, data = d, engine = "BayesX", family = gaussian2)



d <- dgp_beta(
  200,
  mu = list(const = -0.5, type = list(c("pick", "const"))),
  sigma2 = list(const = 0.01, type = list("quadratic", "const"))
)

no <- gaussian2.bamlss()
no$iwls <- NULL

f <- list(
  y ~ s(mu.x11, k = 40),
  sigma ~ s(sigma2.x11)
)

b <- bamlss(f, data = d, family = tF(BE), update = "iwls")


## Tensor tests
test2<-function(u,v,w,sv=0.3,sw=0.4)  
{ ((pi**sv*sw)*(1.2*exp(-(v-0.2)^2/sv^2-(w-0.3)^2/sw^2)+
  0.8*exp(-(v-0.7)^2/sv^2-(w-0.8)^2/sw^2)))*(u-0.5)^2*20
}
n <- 500
v <- runif(n);w<-runif(n);u<-runif(n)
f <- test2(u,v,w)
y <- f + rnorm(n)*0.2
# tensor product of 2D Duchon spline and 1D cr spline
m <- list(c(1,.5),rep(0,0)) ## example of list form of m
b <- gam(y~te(v,w,u,k=c(30,5),d=c(2,1),bs=c("ds","cr"),m=m))
op <- par(mfrow=c(2,2))
vis.gam(b,cond=list(u=0),color="heat",zlim=c(-0.2,3.5))
vis.gam(b,cond=list(u=.33),color="heat",zlim=c(-0.2,3.5))
vis.gam(b,cond=list(u=.67),color="heat",zlim=c(-0.2,3.5))
vis.gam(b,cond=list(u=1),color="heat",zlim=c(-0.2,3.5))
par(op)

b0 <- bamlss(y ~ te(v,w,u,k=c(30,5),d=c(2,1),bs=c("ds","cr"),m=m), method = c("MP2", "MCMC"), propose = "iwls0")

b1 <- bamlss(y ~ te(v,w,u,k=c(30,5),d=c(2,1),bs=c("ds","cr"),m=m), method = "backfitting", update = "iwls")
b1 <- bamlss(y ~ te(v,w,u,k=c(30,5),d=c(2,1),bs=c("ds","cr"),m=m), method = "backfitting", update = "optim2")


b2 <- bamlss(y ~ rs(s(v,w,bs="ds",k=30), s(u,bs="cr",k=5)), method = "backfitting")





library("bamlss")

set.seed(111)
dat <- gamSim(4)
dat <- dat[, c("y", "x0", "fac")]
dat <- cbind(dat, as.data.frame(model.matrix(~ -1 + fac, data = dat)))

b0 <- bamlss(y ~ x0 + x0:fac, data = dat, engine = "BayesX")
b1 <- bamlss(y ~ x0 + x0:fac, data = dat, method = "MP2")
b2 <- gam( y ~ x0 + x0:fac, data = dat)

summary(b0)
summary(b1)
summary(b2)

nd <- data.frame(x0 = rep(0.5, 3), fac = dat$fac[1])

nd$p0 <- predict(b0, newdata = nd, model = "mu")
nd$p1 <- predict(b1, newdata = nd, model = "mu")
nd$p2 <- predict(b2, newdata = nd) 

print(nd)



## Boost
     set.seed(1907)
     x1 <- rnorm(1000)
     x2 <- rnorm(1000)
     x3 <- rnorm(1000)
     x4 <- rnorm(1000)
     x5 <- rnorm(1000)
     x6 <- rnorm(1000)
     mu    <- 1.5 +1 * x1 +0.5 * x2 -0.5 * x3 -1 * x4
     sigma <- exp(-0.4 * x3 -0.2 * x4 +0.2 * x5 +0.4 * x6)
     y <- numeric(1000)
     for( i in 1:1000)
         y[i] <- rnorm(1, mean = mu[i], sd = sigma[i])
     dat <- data.frame(x1, x2, x3, x4, x5, x6, y)
     
     ### linear model with y ~ . for both components: 400 boosting iterations
     model <- glmboostLSS(y ~ ., families = GaussianLSS(), data = dat,
                          control = boost_control(mstop = 400),
                          center = TRUE)

b <- bamlss(y ~ ., ~., data = dat, sampler = NULL, optimizer = boost0)





##################################################
## NLS ###########################################
##################################################
data("O2K", package = "nlstools")

b <- gam(VO2 ~ s(t, k = 10,bs="ps"), data = O2K)

sobj <- smoothCon(s(t, k = 5, bs = "ps"), O2K, NULL, absorb.cons = TRUE)[[1]]
X <- sobj$X
S <- sobj$S[[1]]
k <- ncol(X)

ff <- function(beta, X) {
  f <- (X %*% beta[2:(k + 1)]) / exp(1 + X %*% beta[(k + 2):(2 * k + 1)])
  f <- f - mean(f)
  drop(beta[1] + f)
}

objfun <- function(theta, X, y, S, lambda) {
  sd <- theta[1]
  beta <- theta[-1]
  f <- ff(beta, X)
  bu <- beta[2:(k + 1)]
  bl <- beta[(k + 2):(2 * k + 1)]
  ll <- sum(dnorm(y, mean = f, sd = exp(sd), log = TRUE))
  -1 * ll
}

start <- c(0, 0, rep(0, ncol(X) * 2))

opt <- optim(start, objfun, X = X, y = O2K$VO2, S = S, lambda = c(0.00001, 0.000000000000001), method = "BFGS")

plot(O2K, ylim = range(c(fitted(b), O2K$V02, ff(opt$par[-1], X))))
lines(fitted(b) ~ O2K$t)
lines(ff(opt$par[-1], X) ~ O2K$t, col = "red")

