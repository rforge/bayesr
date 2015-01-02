library("bamlss")

## pick function
f <- simfun(type = "complicated")

set.seed(111)
n <- 1000
dat <- data.frame("x1" = sort(runif(n, 0, 1)))
dat$y <- with(dat, 1.2 + 0.4 * x1 + rnorm(n, sd = exp(-1 + 0.3 * x1)))

theta <- list("mu" = rep(0, 2), "sigma" = rep(0, 2))
d <- list("mu" = cbind(1, dat$x1), "sigma" = cbind(1, dat$x1))

logLik <- function(mu, sigma, ...) {
  dnorm(dat$y, mean = d$mu %*% mu, sd = exp(d$sigma %*% sigma), log = TRUE)
}

b0 <- gmcmc(logLik, theta = theta, propose = gmcmc_propose_default)
b1 <- gmcmc(logLik, theta = theta, propose = gmcmc_slice)
b2 <- gmcmc(logLik, theta = theta, propose = c(gmcmc_slice, gmcmc_propose_default))

apply(b0, 2, mean)
apply(b1, 2, mean)
apply(b2, 2, mean)


library("MCMCpack")

logitfun <- function(beta, y, X){
  eta <- X %*% beta
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


##  negative binomial regression with an improper unform prior
## X and y are passed as args to MCMCmetrop1R
negbinfun <- function(theta, y, X){
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

b <- gmcmc(negbinfun, theta = list(theta = c(0, 0, 0, 0)),
  y = yy, X = XX,
  n.iter = 35000, burnin = 1000, thin = 1)

