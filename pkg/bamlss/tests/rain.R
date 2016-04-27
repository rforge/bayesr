library("bamlss")

if(file.exists("~/data/homstart.rda")) {
  load("~/data/homstart.rda")
} else {
  dpath <- system.file(package = "bamlss", "data")
  if(!file.exists(file <- file.path(dpath, "homstart.rda"))) {
    homstart_data(dir = dirname(file), load = TRUE)
  } else load(file)
  file.copy(file, file.path("~/data/homstart.rda"))
}

homstart$raw[homstart$raw < 0] <- 0

f <- list(
  "mu" = raw ~ s(day,bs="cc"),
  "sigma" = ~ s(day,bs="cc"),
  "alpha" = ~ 1
)

for(i in levels(homstart$id)) {
  td <- subset(homstart, id == i)

  b <- bamlss(f, data = td, family = "pcnorm",
    binning = TRUE, gam.side = FALSE, before = TRUE)

  plot(b)
}


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
f <- list(
  rain ~ s(sqrtensmean),
  sigma ~ s(sqrtenssd),
  alpha ~ 1
)
b1 <- bamlss(f, data = RainIbk, family = "pcnorm", sampler = FALSE)

f <- list(
  I(rain^(1/2)) ~ s(sqrtensmean),
  sigma ~ s(sqrtenssd),
  alpha ~ s(sqrtensmean) + s(sqrtenssd)
)
b2 <- bamlss(f[1:2], data = RainIbk, family = gF(cnorm), no.bfit = TRUE, no.mcmc = TRUE)

f1 <- gF2(pcnorm)
f2 <- gF2(cnorm)

sum(f1$d(RainIbk$rain, f1$map2par(fitted(b1)), log = TRUE))
a <- ifelse(RainIbk$rain <= 0, 0 , - log(2) + (1 / 2 - 1) * log(RainIbk$rain))
sum(f2$d(RainIbk$rain^(1/2), f2$map2par(fitted(b2)), log = TRUE) + a)

plot(c(b1, b2))


RainIbk$p <- exp(predict(b1, model = "alpha"))

b2 <- bamlss(I(rain^(1/p)) ~ s(sqrtensmean), ~ s(sqrtenssd), data = RainIbk, family = cnorm,
  n.iter = 1200, burnin = 200, thin = 10)
b3 <- bamlss(I(rain^(1/2)) ~ s(sqrtensmean), ~ s(sqrtenssd), data = RainIbk, family = cnorm,
  n.iter = 1200, burnin = 200, thin = 1)

f <- gF(cnorm)
bf <- fitted(b2)
sum(with(RainIbk, ifelse(rain <= 0,
  pnorm(0, bf$mu, bf$sigma),
  dnorm(rain^(1/p), bf$mu, bf$sigma, log = TRUE) - log(p) + (1/p-1) * log(rain))))
bf <- fitted(b3)
sum(with(RainIbk, ifelse(rain <= 0,
  pnorm(0, bf$mu, bf$sigma),
  dnorm(rain^(1/p), bf$mu, bf$sigma, log = TRUE) - log(p) + (1/p-1) * log(rain))))

par(mfrow = c(1, 2))
with(subset(RainIbk, rain > 0), hist(rain^(1/p), freq = FALSE))
with(subset(RainIbk, rain > 0), hist(rain^(1/2), freq = FALSE))


boxcox <- function(x, lambda = 1)
{
  i <- x > 0
  x[i] <- if(lambda != 0) (x[i]^lambda - 1) / lambda else log(x[i])
  x
}

alpha <- predict(b2, model = "alpha", type = "parameter")

b1 <- bamlss(sqrt(rain) ~ s(sqrtensmean), ~ s(sqrtenssd), data = RainIbk, family = cnorm)

RainIbk$f_crch <- predict(b0)
RainIbk$f_bamlss <- predict(b1, model = "mu", intercept = TRUE)

plot(sqrt(rain) ~ sqrtensmean, data = RainIbk, pch = ".") 
plot2d(f_crch + f_bamlss ~ sqrtensmean, data = RainIbk, add = TRUE,
  col.lines = c("red", "green"))
legend("topleft", c("crch", "bamlss"), lwd = 1, col = c("red", "green"))



load("~/svn/SnowSafeFX/data/ehyddata.rda")

dat$obs[dat$obs < 0] <- 0
dat$yday <- as.POSIXlt(dat$date)$yday
dat$year <- as.POSIXlt(dat$date)$year + 1900
dat <- subset(dat, year > 2010 & year <= 2011)

f <- list(
  sqrt(obs) ~ te(yday,lon,lat, bs=c("cc","tp"), d=c(1, 2)) + s(alt) + s(lon, lat),
            ~ s(yday, bs = "cc") + s(alt) + s(lon, lat)
)

b1 <- bamlss(f, data = dat, family = "cnorm", n.iter = 120, burnin = 20, thin = 1, binning = TRUE, before = FALSE)


set.seed(111)

n <- 100
d <- data.frame("x" = runif(n))
d$y <- -0.2 + 0.4 * d$x + rnorm(n, sd = 0.3)
d$yobs <- ifelse(d$y <= 0, 0, d$y)

b <- bamlss0(yobs ~ x, data = d, family = gF(cens, left = 0))

coef(b)
plot(b, which = 3:4)

set.seed(111)

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
d <- subset(d, rain > 0)

b <- bamlss(rain ~ s2(lon,lat,year,time,bs="str"), data = d)



f <- list(
  rain ~ ti(time) + ti(time,lon,lat,bs=c("cc","tp"),d=c(1,2), k = c(5, 10)),
  ~ 1
)

b <- bamlss(f, data = d, n.iter = 1200, burnin = 200, thin = 1,
  binning = TRUE, before = TRUE, gam.side = FALSE, family = "cnorm")


