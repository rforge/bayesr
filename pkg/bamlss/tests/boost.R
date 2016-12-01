library("bamlss")

n <- 1000
d <- data.frame(
  "x1" = runif(n),
  "x2" = runif(n),
  "x3" = runif(n),
  "x4" = runif(n),
  "x5" = runif(n),
  "x6" = runif(n)
)
d$eta_mu <- with(d, 1.2 + x1 - 0.5 * x2 + sin(scale2(x3, -3, 3)))
d$eta_sigma <- with(d, 0.5 * x1 + cos(scale2(x6, -3, 3)))
d$y <- rnorm(n, mean = d$eta_mu, sd = exp(d$eta_sigma))

f <- list(
  y ~ x1 + x2 + x3 + x4 + x5 + x6 + s(x1) + s(x2) + s(x3) + s(x4) + s(x5) + s(x6),
  sigma ~ x1 + x2 + x3 + x4 + x5 + x6 + s(x1) + s(x2) + s(x3) + s(x4) + s(x5) + s(x6)
)

b <- bamlss(f, data = d, optimizer = boost, sampler = FALSE, scale.d = TRUE)

load("cars.rda")
f <- price ~ poly(age, 3) + poly(kilometer, 3) +
  poly(TIA, 3) + abs + sunroof - 1
f <- update.formula(f, ~ .^2)

b <- bamlss(f, data = cars, optimizer = boost, sampler = FALSE, scale.d = TRUE)


d <- list()
n <- 100
for(i in 1:n)
  d[[paste("x", i, sep = "")]] <- runif(500)
d <- as.data.frame(d)

j <- sample(1:ncol(d), size = 8, replace = FALSE)
coef <- rep(0, ncol(d))
coef[j] <- runif(length(j), -1, 1)
d$y <- as.matrix(d) %*% coef + rnorm(500, sd = 0.3)

f <- paste(names(d)[-ncol(d)], collapse = "+")
f <- as.formula(paste("y~", f, sep = ""))

b <- bamlss(f, data = d, optimizer = boost, sampler = FALSE, scale.d = TRUE)


d <- GAMart()
b <- bamlss(num ~ x1 + x2 + x3 + s(x1) + s(x2) + s(x3), data = d, optimizer = boost, sampler = FALSE, scale.d = TRUE)


## Rain example.
homstart_data()
homstart$raw[homstart$raw < 0] <- 0

f <- list(
  raw ~ s(day,k=20),
  sigma ~ s(day,k=20),
  lambda ~ s(day,k=20)
)

ff <- bamlss:::pcnorm_bamlss()

nd <- data.frame("day" = 1:365)
for(i in levels(homstart$id)) {
  cat("Station", i, "\n")
  d <- subset(homstart, id == i)
  if(nrow(d) < 100) next
  b <- bamlss(f, data = d, family = ff, optimizer = boost, sampler = FALSE, maxit = 3000)
  p <- predict(b, newdata = nd)
  names(p) <- paste("s", i, ".", names(p), sep = "")
  nd <- cbind(nd, do.call("cbind", p))
}


## Lasso.
n <- 1000
d <- data.frame(
  "x1" = runif(n, -1, 1),
  "x2" = runif(n, -1, 1),
  "x3" = runif(n, -1, 1),
  "x4" = runif(n, -1, 1),
  "x5" = runif(n, -1, 1),
  "x6" = runif(n, -1, 1)
)
d$eta_mu <- with(d, 1.2 + x1 - 0.5 * x2 + x6)
d$eta_sigma <- with(d, 0.5 * x1  - 0.8 * x4)
d$y <- rnorm(n, mean = d$eta_mu, sd = exp(d$eta_sigma))

f <- list(
  y ~ x1 + x2 + x3 + x4 + x5 + x6,
  sigma ~ x1 + x2 + x3 + x4 + x5 + x6
)

b <- bamlss(f, data = d, lasso = TRUE, propose = "iwls", n.iter = 12000, burnin = 2000, thin = 10)


dgp <- function(n = 100, p = 10,
  sd = 0.3, sigma = diag(p), prob = c(0.7, 0.3), kmax = NULL)
{
  require("mvtnorm")

  beta <- if(is.null(kmax)) {
    sample(0:1, size = p, replace = TRUE, prob = prob)
  } else {
    sample(c(rep(0, p - kmax), rep(1, kmax)))
  }

  X <- rmvnorm(n, rep(0, p), sigma = sigma)
  colnames(X) <- paste("x", 1:ncol(X), sep = "")
  f <- X %*% beta
  y <- f + rnorm(n, sd = sd)

  d <- data.frame("y" = y, "f" = f)
  d <- cbind(d, as.data.frame(X))

  attr(d, "beta") <- beta
  return(d)
}

d <- dgp(n = 300, p = 20, sd = 0.1)

f <- grep("x", colnames(d), value = TRUE)
f <- paste(f, collapse = "+")
f <- as.formula(paste("", f, sep = "~"))

b <- bamlss(y ~ la(f), data = d, sampler = FALSE)


b0 <- bamlss(f, data = d, sampler = FALSE)
b1 <- bamlss(f, data = d, lasso = TRUE, sampler = FALSE, criterion = "BIC", sep.lasso = FALSE)
b2 <- bamlss(f, data = d, lasso = TRUE, sampler = FALSE, criterion = "BIC", sep.lasso = TRUE)

cb0 <- coef(b0, hyper = FALSE)
cb0 <- round(cb0[!grepl("Intercept", names(cb0))], 4)

cb1 <- coef(b1, hyper = FALSE)
cb1 <- round(cb1[!grepl("Intercept", names(cb1))], 4)

cb2 <- coef(b2, hyper = FALSE)
cb2 <- round(cb2[!grepl("Intercept", names(cb2))], 4)

cbind("true" = attr(d, "beta"), "fixed" = cb0 != 0, "lasso" = cb1 != 0, "sep.lasso" = cb2 != 0)



data(rent99, package = "gamlss.data")
rent99$district <- as.factor(rent99$district)

f <- rentsqm ~ la(location) + la(yearc) + la(area) +
  la(bath) + la(kitchen) + la(cheating) + la(district)

b <- bamlss(f, data = rent99, scale.d = TRUE, criterion = "BIC")

