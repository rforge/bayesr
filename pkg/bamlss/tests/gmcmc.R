library("bamlss")

## pick function
f <- simfun(type = "pick")

set.seed(111)
n <- 1000
dat <- data.frame("x1" = sort(runif(n, 0, 1)))
dat$y <- with(dat, 1.2 + f(x1) + rnorm(n, sd = exp(-1 + 0.3 * x1)))

theta <- list("mu" = c(-10, 20), "sigma" = rep(0, 2))
d <- list("mu" = cbind(1, dat$x1), "sigma" = cbind(1, dat$x1))

## b <- bamlss0(y ~ s(x1), ~ s(x1), data = dat, n.iter = 1200, burnin = 200, thin = 1)

logLik <- function(mu, sigma, ...) {
  mu <- unlist(mu)
  sigma <- unlist(sigma)
  ll <- dnorm(dat$y, mean = d$mu %*% mu, sd = exp(d$sigma %*% sigma), log = TRUE)
}

b0 <- gmcmc(logLik, theta = theta, propose = gmcmc_mvnorm)
b1 <- gmcmc(logLik, theta = theta, propose = gmcmc_slice)
b2 <- gmcmc(logLik, theta = theta, propose = c(gmcmc_slice, gmcmc_mvnorm))
b3 <- gmcmc(logLik, theta = theta, propose = gmcmc_newton)

g <- list(
  "mu" = function(mu, sigma) {
    mu <- d$mu %*% unlist(mu)
    sigma <- exp(d$sigma %*% unlist(sigma))
    score <- snorm(dat$y, mean = mu, sd = sigma, which = "mu")
    drop(t(d$mu) %*% score)
  },
  "sigma" = function(mu, sigma) {
    mu <- d$mu %*% unlist(mu)
    sigma <- exp(d$sigma %*% unlist(sigma))
    score <- snorm(dat$y, mean = mu, sd = sigma, which = "sigma")
    score <- score * sigma
    drop(t(d$sigma) %*% score)
  }
)

b4 <- gmcmc(logLik, theta = theta, propose = gmcmc_newton, gradient = g)

rbind(
  apply(b0, 2, mean),
  apply(b1, 2, mean),
  apply(b2, 2, mean),
  apply(b3, 2, mean),
  apply(b4, 2, mean)
)


library("MCMCpack")

logitfun <- function(beta, y, X){
  eta <- X %*% unlist(beta)
  p <- 1.0/(1.0+exp(-eta))
  sum( y * log(p) + (1-y)*log(1-p) )
}
        
x1 <- rnorm(1000)
x2 <- rnorm(1000)
Xdata <- cbind(1,x1,x2)
p <- exp(.5 - x1 + x2)/(1+exp(.5 - x1 + x2))
yvector <- rbinom(1000, 1, p)

post.samp <- MCMCmetrop1R(logitfun, theta.init = c(0,0,0),
  X = Xdata, y = yvector,
  thin = 1, mcmc = 40000, burnin = 500,
  tune = c(1.5, 1.5, 1.5),
  verbose = 500, logfun = TRUE)

b <- gmcmc(logitfun, theta = list(beta = c(0, 0, 0)),
  y = yvector, X = Xdata,
  n.iter = 12000, burnin = 2000, thin = 10,
  propose = gmcmc_slice)


## Negative binomial regression with an improper unform prior
## X and y are passed as args to MCMCmetrop1R.
negbinfun <- function(theta, y, X) {
  theta <- unlist(theta)
  k <- length(theta)
  beta <- theta[1:(k - 1)]
  alpha <- exp(theta[k])
  mu <- exp(X %*% beta)
  log.like <- sum(
    lgamma(y + alpha) - lfactorial(y) - lgamma(alpha) +
    alpha * log(alpha / (alpha + mu)) +
    y * log(mu / (alpha + mu))
  )
}
     
n <- 1000
x1 <- rnorm(n)
x2 <- rnorm(n)
XX <- cbind(1, x1, x2)
mu <- exp(1.5 + x1 + 2*x2) * rgamma(n, 1)
yy <- rpois(n, mu)
     
post.samp <- MCMCmetrop1R(negbinfun, theta.init = c(0, 0, 0, 0), y = yy, X = XX,
  thin = 1, mcmc = 35000, burnin = 1000,
  tune = 1.5, verbose = 500, logfun = TRUE,
  seed = list(NA, 1))

b <- gmcmc(negbinfun, theta = c(0, 0, 0, 0),
  y = yy, X = XX,
  n.iter = 35000, burnin = 6000, thin = 1)

logfun <- function(x) {
  dnorm(unlist(x), mean = 5, sd = 0.1, log = TRUE)
}
b <- gmcmc(logfun, theta = list(x = 5.3), n.iter = 1200, burnin = 200, thin = 1)


## Example taken from
## http://www.dme.ufrj.br/mcmc/Example73-R.txt
loglik = function(theta) {
  theta <- unlist(theta)
  mu <- theta[1] / (1 + exp(theta[2]) * (exp(theta[3]) / (1 + exp(theta[3])))^t)
  ll <- sum(dnorm(y, mean = mu, sd = exp(theta[4]), log = TRUE))
}

n = 15
t = 1:n
y = c(16.08, 33.83, 65.80, 97.20, 191.55, 326.20, 386.87, 520.53, 590.03,
  651.92, 724.93, 699.56, 689.96, 637.56, 717.41)

theta <- c(700, 36, 0.5946)
theta <- c(theta, sd(y - theta[1] / (1 + theta[2] * theta[3]^t)))

b <- gmcmc(loglik, theta = theta, propose = gmcmc_newton, n.iter = 12000, burnin = 2000, thin = 10)
apply(b, 2, mean)
theta

fit <- apply(b[, 1:4], 1, function(theta) {
  theta[1] / (1 + exp(theta[2]) * (exp(theta[3]) / (1 + exp(theta[3])))^t)
})

matplot(t, fit, type = "l", col = rgb(0.1, 0.1, 0.1, alpha = 0.01))
points(t, y, col = "blue", lwd = 2)

b1 <- bamlss(I(y / 100) ~ -1 + s2(t, bs = "gc"), ~ 1 + t, method = "MCMC")
b2 <- bamlss(I(y / 100) ~ -1 + s2(t, bs = "gc"), ~ 1 + t, engine = "JAGS",
  n.iter = 12000, burnin = 2000, thin = 10)
b3 <- bamlss0(I(y / 100) ~ s(t), ~ 1 + t)

nd <- data.frame("y" = y / 100, "t" = t)
nd$p1 <- predict(b1, model = 1, FUN = function(x) { quantile(x, probs = c(0.025, 0.975)) })
nd$p2 <- predict(b2, model = 1, FUN = function(x) { quantile(x, probs = c(0.025, 0.975)) })

plot(y ~ t, data = nd)
plot2d(as.matrix(p1) ~ t, data = nd, add = TRUE)
plot2d(as.matrix(p2) ~ t, data = nd, add = TRUE, col.lines = 2)

## Spline example.
library("MASS")
data("mcycle")

S <- smooth.construct(s(times, k = 10), mcycle, NULL)
X <- S$X
K <- S$S[[1]]
a <- b <- 0.0001

logpost = function(mu, sigma) {
  tau21 <- exp(mu$tau2)
  tau22 <- exp(sigma$tau2)

  f1 <- X %*% mu$gamma
  f2 <- exp(X %*% sigma$gamma)

  ll <- sum(dnorm(mcycle$accel, mean = f1, sd = f2, log = TRUE))
  lp1 <- -log(tau21) * S$rank / 2 + drop(-0.5 / tau21 * crossprod(mu$gamma, K) %*% mu$gamma) +
            log((b^a)) - log(gamma(a)) + (-a - 1) * log(tau21) - b / tau21
  lp2 <- -log(tau22) * S$rank / 2 + drop(-0.5 / tau22 * crossprod(sigma$gamma, K) %*% sigma$gamma) +
            log((b^a)) - log(gamma(a)) + (-a - 1) * log(tau22) - b / tau22

  ll + lp1 + lp2
}

theta <- list(
  mu = list("gamma" = rep(0, ncol(X)), "tau2" = 0),
  sigma = list("gamma" = rep(0, ncol(X)), "tau2" = 0)
)

b <- gmcmc(logpost, theta = theta, propose = gmcmc_newton, n.iter = 1200, burnin = 200, thin = 10,
  adapt = 5)

nb <- colnames(b)
fit1 <- apply(b[, grep("mu.gamma[", nb, fixed = TRUE)], 1, function(g) {
  X %*% g
})
fit2 <- apply(b[, grep("sigma.gamma[", nb, fixed = TRUE)], 1, function(g) {
  X %*% g
})

par(mfrow = c(2, 1))
plot(mcycle, col = "blue", lwd = 2)
matplot(mcycle$times, fit1, type = "l", col = rgb(0.1, 0.1, 0.1, alpha = 0.01), add = TRUE)
f <- t(apply(fit1, 1, quantile, probs = c(0.05, 0.5, 0.95)))
matplot(mcycle$times, f, type = "l", lty = c(2, 1, 2), col = 1, add = TRUE)
matplot(mcycle$times, fit2, type = "l", col = rgb(0.1, 0.1, 0.1, alpha = 0.01))
f <- t(apply(fit2, 1, quantile, probs = c(0.05, 0.5, 0.95)))
matplot(mcycle$times, f, type = "l", lty = c(2, 1, 2), col = 1, add = TRUE)


## Simple example.
y <- rnorm(500, mean = 10, sd = 2)

logpost = function(theta) {
  theta <- unlist(theta)
  ll <- sum(dnorm(y, mean = theta[1], sd = exp(theta[2]), log = TRUE))
}

g <- list("theta" = function(theta) {
  theta <- unlist(theta)
  mu <- theta[1]
  sigma <- exp(theta[2])
  score <- snorm(y, mean = mu, sd = sigma)
  score[, "dsigma"] <- score[, "dsigma"] * sigma
  score <- colSums(score)
  return(score)
})

h <- list("theta" = function(theta) {
  theta <- unlist(theta)
  mu <- theta[1]
  sigma <- exp(theta[2])
  score <- snorm(y, mean = mu, sd = sigma, which = "sigma")
  hess <- hnorm(y, mean = mu, sd = sigma,
    which = c("mu", "mu.sigma", "sigma.mu", "sigma"))
  hess[, "dmu.dsigma"] <- hess[, "dsigma.dmu"] * sigma
  hess[, "d2sigma"] <- hess[, "d2sigma"] * sigma^2 + sigma * score
  hess <- colSums(hess)
  hess <- matrix(hess, 2, 2, byrow = TRUE)
  rownames(hess) <- colnames(hess) <- c("mu", "sigma")
  return(-1 * hess)
})

f <- list("theta" = function(theta) {
  theta <- unlist(theta)
  mu <- theta[1]
  sigma <- exp(theta[2])
  fish <- fnorm(y, mean = mu, sd = sigma,
    which = c("mu", "mu.sigma", "sigma.mu", "sigma"))
  fish[, "d2sigma"] <- fish[, "d2sigma"] * sigma^2
  fish <- colSums(fish)
  fish <- matrix(fish, 2, 2, byrow = TRUE)
  rownames(fish) <- colnames(fish) <- c("mu", "sigma")
  return(fish)
})

theta <- c("mu" = 10, "sigma" = 0.6)

b <- gmcmc(logpost, theta = theta, propose = gmcmc_newton, gradient = g, hessian = h)

apply(b, 2, mean)

