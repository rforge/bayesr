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
f <- price ~ -1 + poly(age, 3) + poly(kilometer, 3) +
  poly(TIA, 3) + abs + sunroof
f <- update.formula(f, ~ .^2)

b <- bamlss(f, data = cars, optimizer = boost, sampler = FALSE, scale.d = TRUE)
