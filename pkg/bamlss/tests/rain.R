library("bamlss")

if(file.exists("~/tmp/homstart.rda")) {
  load("~/tmp/homstart.rda")
} else {
  dpath <- system.file(package = "bamlss", "data")
  if(!file.exists(file <- file.path(dpath, "homstart.rda"))) {
    homstart_data(dir = dirname(file), load = TRUE)
  } else load(file)
  file.copy(file, file.path("~/tmp/homstart.rda"))
}

homstart$raw[homstart$raw == 0] <- 0.05
homstart$raw[homstart$raw < 0] <- 0

rain <- subset(homstart, year >= 2008)
rain2 <- subset(homstart, year >= 2008 & raw > 0)

f <- list(
  sqrt(raw) ~ s(day, bs = "cc") + s(elevation) + s(long, lat),
  ~ s(day, bs = "cc") + s(elevation) + s(long, lat)
)

b1 <- bamlss(f, data = rain2, family = gF(cens, left = 0),
  method = c("backfitting", "MCMC"), update = "iwls", propose = "iwls",
  n.iter = 5000, burnin = 1000, thin = 10, cores = 3)

library("truncreg")

set.seed(1071)
n <- 10000
sigma <- 4
alpha <- 2
beta <- 1
x <- rnorm(n, mean = 0, sd = 2)
eps <- rnorm(n, sd = sigma)
y <- alpha + beta * x + eps
d <- data.frame(y = y, x = x)

## truncated response
d$yt <- ifelse(d$y > 1, d$y, NA)

b0 <- truncreg(yt ~ x, data = d, point = 1, direction = "left")

b <- bamlss(yt ~ x, data = d, family = gF(trunc, point = 1))
mf <- model.frame(b)
eta <- fitted(b)
dy <- gF(truncreg, point = 1)$d(mf$yt, eta)
hist(mf$yt, freq = FALSE)
lines(dy[order(mf$yt)] ~ mf$yt[order(mf$yt)], col = 2)


## CRCH tests
library("crch")
data("RainIbk")

## mean and standard deviation of square root transformed ensemble forecasts
RainIbk$sqrtensmean <- apply(sqrt(RainIbk[,grep('^rainfc',names(RainIbk))]), 1, mean)
RainIbk$sqrtenssd <- apply(sqrt(RainIbk[,grep('^rainfc',names(RainIbk))]), 1, sd)
   
## left censored regression model with censoring point 0 and 
## conditional heteroscedasticy:
b0 <- crch(sqrt(rain) ~ sqrtensmean|sqrtenssd, data = RainIbk, dist = "gaussian",  left = 0)

## now with bamlss
b1 <- bamlss(sqrt(rain) ~ sqrtensmean, ~ sqrtenssd, data = RainIbk, family = cens,
  method = c("backfitting", "MCMC"), update = "iwls", propose = "iwls",
  n.iter = 1200, burnin = 200, thin = 2)

RainIbk$f_crch <- predict(b0)
RainIbk$f_bamlss <- predict(b1, model = "mu", intercept = TRUE)

plot(sqrt(rain) ~ sqrtensmean, data = RainIbk, pch = ".") 
plot2d(f_crch + f_bamlss ~ sqrtensmean, data = RainIbk, add = TRUE,
  col.lines = c("red", "green"), lwd = c(6, 3))
legend("topleft", c("crch", "bamlss"), lwd = 1, col = c("red", "green"))

