library("bamlss")
library("maptools")
library("rgdal")
library("raster")

if(file.exists("homstart.rda")) {
  load("homstart.rda")
} else {
  dpath <- system.file(package = "bamlss", "data")
  if(!file.exists(file <- file.path(dpath, "homstart.rda"))) {
    homstart_data(dir = dirname(file), load = TRUE)
  } else load(file)
  file.copy(file, file.path("homstart.rda"))
}

homstart$raw[homstart$raw < 0] <- 0

if(!file.exists("rainmodel.rda")) {
  homstart$raw[homstart$raw < 0] <- 0
  homstart2 <- subset(homstart, year >= 1979)

  f <- list(
    "mu" = sqrt(raw) ~ elevation + ti(day,bs="cc",k=20) + ti(long,lat,bs="tp",d=2,k=30) +
      ti(day,long,lat,bs=c("cc","tp"),d=c(1,2),mp=FALSE),
    "sigma" = ~ elevation + ti(day,bs="cc",k=20) + ti(long,lat,bs="tp",d=2,k=30) +
      ti(day,long,lat,bs=c("cc","tp"),d=c(1,2),mp=FALSE)
  )

  rainmodel <- bamlss(f, data = homstart2, family = "cnorm",
    binning = TRUE, before = TRUE, gam.side = FALSE,
    samplestats = FALSE, results = FALSE,
    eps = 0.001, n.iter = 6000, burnin = 4000, thin = 10, cores = 7)

  save(rainmodel, file = "rainmodel.rda")
}

if(FALSE) {
load("rainmodel.rda")

expCens <- function(mu, sigma) {
  pnorm(mu / sigma) * (mu + sigma * dnorm(mu / sigma) / pnorm(mu / sigma))
}

expSample <- function(mu, sigma) {
  n <- 1000
  N <- length(mu)
  rval <- sapply(1:N, function(i) {
    rs <- rnorm(n, mean = mu[i], sd = sigma[i])
    rs[rs < 0] <- 0
    rs <- mean(rs^2)
    rs
  })
  unlist(rval)
}

load("data/AUT_etopo.rda")
load("data/colors.rain.rda")

austria <- readShapePoly(file.path("data", "shp", "at"),
  proj4string = CRS("+proj=longlat +datum=WGS84"))
austria <- as(austria, "SpatialPolygons")
austria2 <- unionSpatialPolygons(austria, as.factor(rep(1, nrow(coordinates(austria)))))
co <- do.call("rbind", sp2bnd(austria2))
n <- 1000
nd2 <- expand.grid("long" = seq(min(co[, 1]), max(co[, 2]), length = n),
  "lat" = seq(min(co[, 1]), max(co[, 2]), length = n))
i <- point.in.polygon(nd2$long, nd2$lat, co[, 1], co[, 2])
nd2 <- subset(nd2, i > 0)
nd2$fmu <- predict(b, newdata = nd2, model = "mu", term = "s(long,lat)", intercept = FALSE)
nd2$fsigma <- predict(b, newdata = nd2, model = "sigma", term = "s(long,lat)", intercept = FALSE)
nd2$elevation <- extract(etopo, cbind(nd2$long, nd2$lat))

nd2$day <- 10
nd2$fmu2 <- predict(b, newdata = nd2, model = "mu", intercept = TRUE)
nd2$fsigma2 <- predict(b, newdata = nd2, model = "sigma", intercept = TRUE)
nd2$rain10 <- expSample(nd2$fmu2, nd2$fsigma2)

nd2$day <- 192
nd2$fmu2 <- predict(b, newdata = nd2, model = "mu", intercept = TRUE)
nd2$fsigma2 <- predict(b, newdata = nd2, model = "sigma", intercept = TRUE)
nd2$rain192 <- expSample(nd2$fmu2, nd2$fsigma2)

nd <- data.frame("day" = sort(unique(homstart$day)))
co <- unique(homstart[, c("long", "lat", "elevation", "id")])
nd3 <- list()
for(i in 1:nrow(co)) {
  cat("Station", i, "\n")
  nd$long <- co[i, "long"]
  nd$lat <- co[i, "lat"]
  nd$elevation <- mean(co$elevation)
  nd[[paste("fmu", i, sep = "")]] <- predict(b, model = "mu", term = "day",
    newdata = nd, intercept = FALSE)
  nd[[paste("fsigma", i, sep = "")]] <- predict(b, model = "sigma", term = "day",
    newdata = nd, intercept = FALSE)
  
  ttime <- subset(homstart, id == co$id[i])
  ttime$day2 <- as.factor(ttime$day)
  fday <- aggregate(sqrt(ttime$raw), by = list(ttime$day2), FUN = function(x) { mean(x, na.rm = TRUE ) })
  names(fday) <- c("day", "mean")
  fday$day <- as.integer(fday$day)
  fday <- fday[order(fday$day), ]
  nd$day <- fday$day
  foo <- function(x) { x }
  pmu <- predict(b, model = "mu", newdata = nd, intercept = TRUE, FUN = foo)
  psigma <- predict(b, model = "sigma", newdata = nd, intercept = TRUE, trans = exp, FUN = foo)
  pred <- NULL
  for(j in 1:ncol(pmu))
    pred <- cbind(pred, expCens(pmu[, j], psigma[, j]))
  pred <- apply(pred, 1, mean, na.rm = TRUE)
  fday$pred <- pred
  attr(fday, "co") <- as.list(co[i, c("long", "lat")])
  nd3[[paste("s", i, sep = "")]] <- fday
}


## Largish version.
png("effects.png", units = "in", res = 150, width = 5.45, height = 3.7 * 1.8)
par(mfrow = c(3, 2), mar = c(1.4, 0.4, 3.1, 0.4), oma = c(4, 4, 0, 2))
redblue1 <- hcl(c(260, 0, 0), c(80, 0, 80), 30, alpha = c(0.15, 0.1, 0.15))
redblue2 <- hcl(c(260, 0, 0), c(80, 0, 80), 30, alpha = c(0.15, 0.1, 0.15) * 3)
latbreaks <- c(46.4, 46.7, 48.2, 48.8)
latcat <- cut(co$lat, breaks = latbreaks)
ylim.mu <- max(abs(range(c(nd[, grep("fmu", names(nd))], nd2$fmu))))
ylim.sigma <- max(abs(range(c(nd[, grep("fsigma", names(nd))], nd2$fsigma))))
ylim.mu <- c(-1, 1) * ylim.mu
ylim.sigma <- c(-1, 1) * ylim.sigma
matplot(nd$day, nd[, grep("fmu", names(nd))], type = "l", lty = 1,
  col = redblue1[latcat], lwd = 2, axes = FALSE, xlab = "", ylab = "",
  ylim = ylim.mu)
legend(80, -0.17, c("Latitude", paste(latbreaks[3:1], latbreaks[4:2], sep = "-")),
  col = c(NA, rev(redblue2)), lty = 1, lwd = 2, bty = "n", cex = 0.8)
abline(h = 0, lty = 2)
Axis(as.Date(c("1970-01-01", "1970-12-31")), side = 1)
axis(2)
box()
mtext("Seasonal effect", side = 2, line = 2.5)
mtext("Time", side = 1, line = 2.5)
mtext(expression(paste("Effects on ", mu)), side = 3, line = 1)
matplot(nd$day, nd[, grep("fsigma", names(nd))], type = "l", lty = 1,
  col = redblue1[latcat], lwd = 2, axes = FALSE, xlab = "", ylab = "",
  ylim = ylim.sigma)
abline(h = 0, lty = 2)
Axis(as.Date(c("1970-01-01", "1970-12-31")), side = 1)
axis(4)
box()
mtext(expression(paste("Effects on ", log(sigma))), side = 3, line = 1)
mtext("Time", side = 1, line = 2.5)

plot(austria2)
xymap(long, lat, fmu, data = nd2, pos = "topleft", layout = FALSE, map = FALSE,
  add = TRUE, color = diverge_hcl, shift = c(0.09, 0.04), distance.labels = 0,
  width = 0.3, symmetric = TRUE, swap = FALSE, range = ylim.mu, digits = 1)
plot(austria, add = TRUE)
box()
axis(1)
axis(2)
mtext("Latitude [deg]", side = 2, line = 3)
legend("bottomleft", "Spatial effect", box.col = NA, bg = NA, cex = 0.9)

plot(austria2)
xymap(long, lat, fsigma, data = nd2, pos = "topleft", layout = FALSE, map = FALSE,
  add = TRUE, color = diverge_hcl, shift = c(0.09, 0.04), distance.labels = 0,
  width = 0.3, symmetric = TRUE, swap = FALSE, range = ylim.sigma, digits = 1)
plot(austria, add = TRUE)
box()
axis(1)
axis(4)
box()

plot(austria2)
xymap(long, lat, rain10, data = nd2, pos = "topleft", layout = FALSE, map = FALSE,
  add = TRUE, color = colors.rain, shift = c(0.09, 0.04), distance.labels = 0,
  width = 0.3, symmetric = FALSE, swap = FALSE, digits = 1, range = c(0, 0.7),
  lrange = round(c(0, max(nd2$rain10)), 1))
plot(austria, add = TRUE)
legend("bottomleft", "Precipitation in mm", box.col = NA, bg = NA, cex = 0.9)
box()
axis(1)
axis(2)
mtext("Longitude [deg]", side = 1, line = 3)
mtext("Latitude [deg]", side = 2, line = 3)
mtext("January", side = 3, line = 0.5)
box()

plot(austria2)
xymap(long, lat, rain192, data = nd2, pos = "topleft", layout = FALSE, map = FALSE,
  add = TRUE, color = colors.rain, shift = c(0.09, 0.04), distance.labels = 0,
  width = 0.3, symmetric = FALSE, swap = FALSE, digits = 1, range = c(0, 2.2),
  lrange = round(c(0, max(nd2$rain192)), 1))
plot(austria, add = TRUE)
box()
axis(1)
axis(4)
mtext("Longitude [deg]", side = 1, line = 3)
mtext("July", side = 3, line = 0.5)
box()
dev.off()


## Small version.
png("effects.png", units = "in", res = 250, width = 6.41, height = 3.1)
par(mfrow = c(1, 2), mar = c(0.1, 0, 0.1, 0), oma = c(3.5, 3.5, 1.5, 0.1))
plot(austria2)
xymap(long, lat, rain10, data = nd2, pos = "topleft", layout = FALSE, map = FALSE,
  add = TRUE, color = colors.rain, shift = c(0.09, 0.04), distance.labels = 0,
  cex.labels = 0.8, length.ticks = 0.5,
  width = 0.3, symmetric = FALSE, swap = FALSE, digits = 1, range = c(0, 0.7),
  lrange = round(c(0, max(nd2$rain10)), 1))
plot(austria, add = TRUE)
for(i in 1:length(nd3)) {
  co2 <- attr(nd3[[i]], "co")
  points(co2$long, co2$lat, pch = 6, cex = 0.3)
}
legend(11.9, 49.792, "Precipitation in mm", box.col = NA, bg = NA, cex = 0.9)
legend("bottomleft", "Available stations", box.col = NA, bg = NA, cex = 0.9, pch = 6)
box()
axis(1)
axis(2)
mtext("Latitude [deg]", side = 2, line = 2.5)
mtext("January", side = 3, line = 0.5)
box()

plot(austria2)
xymap(long, lat, rain192, data = nd2, pos = "topleft", layout = FALSE, map = FALSE,
  add = TRUE, color = colors.rain, shift = c(0.09, 0.04), distance.labels = 0,
  cex.labels = 0.8, length.ticks = 0.5,
  width = 0.3, symmetric = FALSE, swap = FALSE, digits = 1, range = c(0, 2.2),
  lrange = round(c(0, max(nd2$rain192)), 1))
plot(austria, add = TRUE)
for(i in 1:length(nd3)) {
  co2 <- attr(nd3[[i]], "co")
  points(co2$long, co2$lat, pch = 6, cex = 0.3)
}
box()
axis(1)
mtext("July", side = 3, line = 0.5)
box()
mtext("Longitude [deg]", side = 1, line = 2, outer = TRUE)
dev.off()


## Plot deviance from tru mean obs.
plot.mean.fit <- function(stations = 1) {
  ylim <- NULL
  for(i in stations)
    ylim <- c(ylim , range(nd3[[i]][, 2:3]))
  ylim <- range(ylim)
  par(mfrow = c(2, 1))
  plot(1, 1, type = "n", ylim = ylim, xlim = c(1, 365), axes = FALSE,
    xlab = "Time", ylab = "Raw vs. fitted values")
  Axis(as.Date(c("1970-01-01", "1970-12-31")), side = 1)
  axis(2)
  box()
  for(i in stations) {
    lines(mean ~ day, data = nd3[[i]], col = rgb(1, 0, 0, alpha = 0.2), lwd = 2)
    lines(pred ~ day, data = nd3[[i]], col = rgb(0.1, 0.1, 0.1, alpha = 0.2), lwd = 2)
  }
  plot(Austria)
  for(i in stations) {
    co <- attr(nd3[[i]], "co")
    points(co$long, co$lat, pch = 16, cex = 1.5)
    text(co$long, co$lat, i, pos = 3, cex = 1.5)
  }
}

plot.mean.fit()


## Animation.
if(!file.exists("animation")) dir.create("animation")
for(drop in c("fmu", "fsigma", "fmu2", "fsigma2", "rain10", "rain192"))
  nd2[[drop]] <- NULL
for(j in 1:365) {
  cat("day", j, "\n")

  nd2$day <- j
  fmu3 <- predict(b, newdata = nd2, model = "mu", intercept = TRUE)
  fsigma3 <- predict(b, newdata = nd2, model = "sigma", intercept = TRUE)
  prain <- expSample(fmu3, fsigma3)

  fmu4 <- predict(b, newdata = nd2, model = "mu", term = "day", intercept = FALSE)

  png(paste("animation/effects", j, ".png", sep = ""),
    units = "in", res = 120, width = 4.1, height = 8.1)
  par(mfrow = c(3, 1), mar = c(4.1, 4.1, 3.1, 1.1))
  redblue1 <- hcl(c(260, 0, 0), c(80, 0, 80), 30, alpha = c(0.15, 0.1, 0.15))
  redblue2 <- hcl(c(260, 0, 0), c(80, 0, 80), 30, alpha = c(0.15, 0.1, 0.15) * 3)
  latbreaks <- c(46.4, 46.7, 48.2, 48.8)
  latcat <- cut(co$lat, breaks = latbreaks)
  ylim.mu <- max(abs(range(c(nd[, grep("fmu", names(nd))], nd2$fmu))))
  ylim.sigma <- max(abs(range(c(nd[, grep("fsigma", names(nd))], nd2$fsigma))))
  ylim.mu <- c(-1, 1) * ylim.mu
  ylim.sigma <- c(-1, 1) * ylim.sigma
  matplot(nd$day, nd[, grep("fmu", names(nd))], type = "l", lty = 1,
    col = redblue1[latcat], lwd = 2, axes = FALSE, xlab = "", ylab = "",
    ylim = ylim.mu, main = "Spatially varying seasonal effect (1)")
  abline(h = 0, lty = 2)
  abline(v = j, lwd = 3)
  legend("bottom", c("Latitude", paste(latbreaks[3:1], latbreaks[4:2], sep = "-")),
    col = c(NA, rev(redblue2)), lty = 1, lwd = 2, cex = 0.9, box.col = NA, bg = "white")
  Axis(as.Date(c("1970-01-01", "1970-12-31")), side = 1)
  axis(2)
  box()
  mtext(expression(paste("Effects on ", mu)), side = 2, line = 2.5)
  mtext("Time", side = 1, line = 2.5)

  plot(austria2, main = "Spatially varying seasonal effect (2)")
  xymap(long, lat, fmu4, data = nd2, pos = "topleft",
    layout = FALSE, map = FALSE,
    add = TRUE, color = diverge_hcl, shift = c(0.04, 0.04), distance.labels = 0,
    cex.labels = 0.8, length.ticks = 0.5,
    width = 0.25, symmetric = TRUE, swap = FALSE, digits = 1, range = c(-0.9, 0.9),
    lrange = c(-1, 1))
  plot(austria, add = TRUE)
  for(i in 1:length(nd3)) {
    co2 <- attr(nd3[[i]], "co")
    points(co2$long, co2$lat, pch = 6, cex = 0.6)
  }
  legend(11.3, 49.13, expression(paste("Effects on ", mu)), box.col = NA, bg = NA, cex = 0.9)
  legend("bottomleft", "Available stations", box.col = NA, bg = NA, cex = 0.9, pch = 6)
  box()
  axis(1)
  axis(2)
  mtext("Latitude [deg]", side = 2, line = 2.5)
  mtext("Longitude [deg]", side = 1, line = 2.5)
  box()

  plot(austria2, main = "Predicted average precipitation")
  xymap(long, lat, prain, data = nd2, pos = "topleft",
    layout = FALSE, map = FALSE,
    add = TRUE, color = colors.rain, shift = c(0.04, 0.04), distance.labels = 0,
    cex.labels = 0.8, length.ticks = 0.5,
    width = 0.25, symmetric = FALSE, swap = FALSE, digits = 1, range = c(0, 2.2),
    lrange = round(c(0, 3), 1))
  plot(austria, add = TRUE)
  for(i in 1:length(nd3)) {
    co2 <- attr(nd3[[i]], "co")
    points(co2$long, co2$lat, pch = 6, cex = 0.6)
  }
  legend(11.3, 49.13, "Precipitation in mm", box.col = NA, bg = NA, cex = 0.9)
  legend("bottomleft", "Available stations", box.col = NA, bg = NA, cex = 0.9, pch = 6)
  box()
  axis(1)
  axis(2)
  mtext("Latitude [deg]", side = 2, line = 2.5)
  mtext("Longitude [deg]", side = 1, line = 2.5)
  box()
  dev.off()
}

for(i in 1:365) {
  tf <- tempfile()
  file.copy(paste("effects", i, ".png", sep = ""), tf)
  system(paste("convert -resize 300x", tf, tf))
  file.copy(tf, paste("p2/effects", formatC(i, width = 3, flag = 0), ".png", sep = ""))
  unlink(tf)
}

}
