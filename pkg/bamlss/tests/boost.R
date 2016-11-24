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
  "x1" = runif(n, -1, 1)
  "x2" = runif(n, -1, 1)
  "x3" = runif(n, -1, 1)
  "x4" = runif(n, -1, 1)
  "x5" = runif(n, -1, 1)
  "x6" = runif(n, -1, 1)
)
d$eta_mu <- with(d, 1.2 + x1 - 0.5 * x2 + x6)
d$eta_sigma <- with(d, 0.5 * x1  -0.8 * x4)
d$y <- rnorm(n, mean = d$eta_mu, sd = exp(d$eta_sigma))

f <- list(
  y ~ x1 + x2 + x3 + x4 + x5 + x6,
  sigma ~ x1 + x2 + x3 + x4 + x5 + x6
)

b <- bamlss(f, data = d, sampler = FALSE, lasso = TRUE)

