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

homstart$raw[homstart$raw < 0] <- 0

rain2 <- subset(homstart, year >= 2008)

f <- list(
  sqrt(raw) ~ te(day,long,lat, bs=c("cc","tp"), d=c(1,2)) + s(elevation,k=4) + s(long,lat),
  ~ te(day,long,lat, bs=c("cc","tp"), d=c(1,2)) + s(elevation,k=4) + s(long,lat)
)

b2 <- bamlss0(f, data = rain2, family = gF(cens, left = 0),
  binning = FALSE, sampler = NULL, do.optim = TRUE, before = FALSE, maxit = 2)

b1 <- bamlss(f, data = homstart, family = gF(cens, left = 0),
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
ff <- gF2(trunc, point = 1)
dy <- ff$d(mf$yt, ff$map2par(eta))
hist(mf$yt, freq = FALSE)
lines(dy[order(mf$yt)] ~ mf$yt[order(mf$yt)], col = 2)


## CRCH tests
library("crch")
library("bamlss")

data("RainIbk")

## mean and standard deviation of square root transformed ensemble forecasts
RainIbk$sqrtensmean <- apply(sqrt(RainIbk[,grep('^rainfc',names(RainIbk))]), 1, mean)
RainIbk$sqrtenssd <- apply(sqrt(RainIbk[,grep('^rainfc',names(RainIbk))]), 1, sd)
   
## left censored regression model with censoring point 0 and 
## conditional heteroscedasticy:
b0 <- crch(sqrt(rain) ~ sqrtensmean|sqrtenssd, data = RainIbk, dist = "gaussian",  left = 0)

## now with bamlss
b1 <- bamlss(sqrt(rain) ~ s(sqrtensmean), ~ s(sqrtenssd), data = RainIbk, family = cens, n.samples = 0)

RainIbk$f_crch <- predict(b0)
RainIbk$f_bamlss <- predict(b1, model = "mu", intercept = TRUE)

plot(sqrt(rain) ~ sqrtensmean, data = RainIbk, pch = ".") 
plot2d(f_crch + f_bamlss ~ sqrtensmean, data = RainIbk, add = TRUE,
  col.lines = c("red", "green"), lwd = c(6, 3))
legend("topleft", c("crch", "bamlss"), lwd = 1, col = c("red", "green"))



load("~/svn/SnowSafeFX/data/ehyddata.rda")

dat$obs[dat$obs == 0] <- 0.01
dat$yday <- as.POSIXlt(dat$date)$yday

f <- list(
  sqrt(obs) ~ s(yday, bs = "cc") + s(alt) + s(lon, lat),
            ~ s(yday, bs = "cc") + s(alt) + s(lon, lat)
)

b1 <- bamlss(f, data = dat, family = gF(cens, left = 0),
  method = c("backfitting", "MCMC"), update = "iwls", propose = "iwls",
  n.iter = 100, burnin = 0, thin = 1)


set.seed(111)

n <- 100
d <- data.frame("x" = runif(n))
d$y <- -0.2 + 0.4 * d$x + rnorm(n, sd = 0.3)
d$yobs <- ifelse(d$y <= 0, 0, d$y)

b <- bamlss0(yobs ~ x, data = d, family = gF(cens, left = 0))

coef(b)
plot(b, which = 3:4)

n <- 6
years <- 2
co <- expand.grid("lon" = seq(0.001, 1, length = n), "lat" = seq(0.001, 1, length = n))
d <- NULL
for(j in 1:years) {
  for(i in 1:nrow(co)) {
    d <- rbind(d, data.frame(
      "year" = j,
      "time" = 1:365,
      "lon" = co$lon[i],
      "lat" = co$lat[i],
      "id" = i
    ))
  }
}


d$rain <- with(d, sin(scale2(time, -pi, pi)) * scale2(sin(lon), 0.3, 0.5) * scale2(lat, 0.1, 0.5))
co <- unique(d[, c("lon", "lat")])
plot(rain ~ time, data = d, type = "n")
for(j in unique(d$id)) {
  d2 <- subset(d, id == j)
  lines(rain ~ time, data = d2)
}

d$rain <- d$rain + rnorm(n, sd = 0.3)

b <- bamlss(rain ~ ti(time,bs="cc") + ti(lon,lat) + ti(time,lon,lat,bs=c("cc","tp"),d=c(1,2)),
  data = d, n.iter = 200, burnin = 0, thin = 1, propose = "iwls0")


