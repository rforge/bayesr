if(!file.exists("figures"))
  dir.create("figures")

epng <- function(file, height = 4, width = 5, res = 200, bg = "white") {
  png(file, height = height, width = width,
    units = "in", res = res, bg = bg)
}

main <- function(text, side = 3, line = 0.5, cex = 1.2, ...) {
  mtext(text, side = side, line = line, cex = cex, ...)
}

set.seed(123)


##############################
## (1) Precipitation model. ##
##############################
## Austria: http://www.statistik.at/web_de/klassifikationen/regionale_gliederungen/nuts_einheiten/index.html
## AustriaTopo: https://www.ngdc.noaa.gov/mgg/global/global.html
if(!file.exists("figures/rainmodel-effects-predict.png")) {
  library("bamlss")
  library("spatstat")
  library("sp")
  library("maptools")
  library("raster")
  library("rgeos")

  if(file.exists("homstart.rda")) {
    load("homstart.rda")
  } else {
    homstart_data(load = TRUE)
    save(homstart, file = "homstart.rda")
  }

  homstart$raw[homstart$raw < 0] <- 0
  homstart <- subset(homstart, year >= 1979)

  data("Austria", package = "bamlss")
  load("rainmodel.rda")

  nd <- as.data.frame(coordinates(AustriaTopo))
  names(nd) <- c("lon", "lat")
  nd$elevation <- extract(AustriaTopo, cbind(nd$lon, nd$lat))
  nd <- na.omit(nd)

  nd$fmu <- predict(rainmodel, newdata = nd, model = "mu", term = "(lon,lat)", intercept = FALSE)
  nd$fsigma <- predict(rainmodel, newdata = nd, model = "sigma", term = "(lon,lat)", intercept = FALSE)

  expCens <- function(mu, sigma) {
    pnorm(mu / sigma) * (mu + sigma * dnorm(mu / sigma) / pnorm(mu / sigma))
  }

  expSample <- function(mu, sigma, n = 1000) {
    N <- length(mu)
    rval <- sapply(1:N, function(i) {
      rs <- rnorm(n, mean = mu[i], sd = sigma[i])
      rs[rs < 0] <- 0
      rs <- mean(rs^2)
      rs
    })
    unlist(rval)
  }

  nd$day <- 10
  nd$fmu2 <- predict(rainmodel, newdata = nd, model = "mu", intercept = TRUE)
  nd$fsigma2 <- predict(rainmodel, newdata = nd, model = "sigma", intercept = TRUE)
  nd$rain10 <- expSample(nd$fmu2, nd$fsigma2)

  nd$day <- 192
  nd$fmu2 <- predict(rainmodel, newdata = nd, model = "mu", intercept = TRUE)
  nd$fsigma2 <- predict(rainmodel, newdata = nd, model = "sigma", intercept = TRUE)
  nd$rain192 <- expSample(nd$fmu, nd$fsigma2)

  nd2 <- data.frame("day" = 1:365)
  co <- unique(homstart[, c("lon", "lat", "elevation", "id")])
  nd3 <- list()
  for(i in 1:nrow(co)) {
    cat("Station", i, "\n")
    nd2$lon <- co[i, "lon"]
    nd2$lat <- co[i, "lat"]
    nd2$elevation <- mean(co$elevation)
    nd2[[paste("fmu", i, sep = "")]] <- predict(rainmodel, model = "mu", term = "day",
      newdata = nd2, intercept = FALSE)
    nd2[[paste("fsigma", i, sep = "")]] <- predict(rainmodel, model = "sigma", term = "day",
      newdata = nd2, intercept = FALSE)
  
    ttime <- subset(homstart, id == co$id[i])
    ttime$day2 <- as.factor(ttime$day)
    fday <- aggregate(sqrt(ttime$raw), by = list(ttime$day2), FUN = function(x) { mean(x, na.rm = TRUE ) })
    names(fday) <- c("day", "mean")
    fday$day <- as.integer(fday$day)
    fday <- fday[order(fday$day), ]
    nd2$day <- fday$day
    foo <- function(x) { x }
    pmu <- predict(rainmodel, model = "mu", newdata = nd2, intercept = TRUE, FUN = foo)
    psigma <- predict(rainmodel, model = "sigma", newdata = nd2, intercept = TRUE, FUN = foo)
    pred <- NULL
    for(j in 1:ncol(pmu))
      pred <- cbind(pred, expCens(pmu[, j], exp(psigma[, j])))
    pred <- apply(pred, 1, mean, na.rm = TRUE)
    fday$pred <- pred
    attr(fday, "co") <- as.list(co[i, c("lon", "lat")])
    nd3[[paste("s", i, sep = "")]] <- fday
  }

  nd4 <- data.frame("day" = 1:365)
  nd4$fdaymu <- predict(rainmodel, newdata = nd4, model = "mu", term = "(day)", intercept = FALSE)
  nd4$fdaysigma <- predict(rainmodel, newdata = nd4, model = "sigma", term = "(day)", intercept = FALSE)

  colors.rain <- function (n, h = c(-160, -38), c. = c(10, 80), l = c(86, 39), 
    power = c(2.75806451612903, 1), fixup = FALSE, gamma = NULL, 
    alpha = 1, ...) 
  {
    if (!is.null(gamma)) 
        warning("'gamma' is deprecated and has no effect")
    if (n < 1L) 
        return(character(0L))
    h <- rep(h, length.out = 2L)
    c <- rep(c., length.out = 2L)
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, 0, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L], 
        C = c[2L] - diff(c) * rval^power[1L], H = h[2L] - diff(h) * 
            rval), fixup = fixup, ...)
    if (!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
        alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
            width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
  }

  if(!file.exists("rainmodel0.rda")) {
    f <- list(sqrt(raw) ~ 1, sigma ~ 1)
    rainmodel0 <- bamlss(f, data = homstart, family = "cnorm",
      results = FALSE, samplestats = FALSE,
      binning = TRUE, before = TRUE, gam.side = FALSE)
    save(rainmodel0, file = "rainmodel0.rda")
  } else {
    load("rainmodel0.rda")
  }

  stations <- unique(homstart[, c("lon", "lat")])
  orange <- rgb(242, 146, 0, maxColorValue = 255)

  plotUnionAustria <- function(...) {
    plot(unionSpatialPolygons(Austria,
      rep(1L, length = length(Austria))), ...)
  }

  epng("figures/rainmodel-data-stations.png", width = 4.5, height = 3.5)
  par(mar = c(4.1, 4.1, 1.5, 1.5))
  plot(Austria, col = gray(0.9), lwd = 0.3)
  plotUnionAustria(add = TRUE)
  points(stations[, 1], stations[, 2], pch = 16, cex = 1, col = "blue")
  points(stations[, 1], stations[, 2], cex = 1)
  box()
  axis(1)
  axis(2)
  mtext("Longitude [deg]", side = 1, line = 3)
  mtext("Latitude [deg]", side = 2, line = 3)
  legend("topleft", "Meteorological station", pch = 21,
    col = "black",
    pt.bg = "blue",
    box.col = NA, bg = NA, cex = 0.95)
  dev.off()

  epng("figures/rainmodel-data-hist.png", width = 4.5, height = 3.5)
  par(mar = c(4.1, 4.1, 1.5, 1.5))
  mf <- model.frame(rainmodel0)
  y <- mf[["sqrt(raw)"]]
  hist(y, freq = FALSE, main = "", col = gray(0.9), breaks = 0:16 - 0.5,
    xlab = expression(paste("Daily ", sqrt(observations))), ylab = "Density")
  ff <- family(rainmodel0)
  cr <- coef(rainmodel0)
  eta <- list(mu = cr[1], sigma = cr[2])
  y2 <- seq(min(y), max(y), length = 200)
  dy <- ff$d(y = y2, par = ff$map2par(eta))
  lines(dy ~ y2, col = orange, lwd = 2)
  legend("right", "Fitted censored distribution", lwd = 2, col = orange,
    box.col = NA, bg = NA, cex = 0.95)
  box()
  dev.off()

  epng("figures/rainmodel-effects-season-mu.png", width = 4.5, height = 3.5)
  par(mar = c(4.1, 4.1, 1.5, 1.5))
  redblue1 <- hcl(c(260, 0, 0), c(80, 0, 80), 30, alpha = c(0.15, 0.1, 0.15))
  redblue2 <- hcl(c(260, 0, 0), c(80, 0, 80), 30, alpha = c(0.15, 0.1, 0.15) * 3)
  latbreaks <- c(46.4, 46.7, 48.2, 48.8)
  latcat <- cut(co$lat, breaks = latbreaks)
  ylim.mu <- max(abs(range(c(nd2[, grep("fmu", names(nd2))], nd2$fmu))))
  ylim.sigma <- max(abs(range(c(nd2[, grep("fsigma", names(nd2))], nd2$fsigma))))
  ylim.mu <- c(-1, 1) * ylim.mu
  ylim.sigma <- c(-1, 1) * ylim.sigma
  matplot(nd2$day, nd2[, grep("fmu", names(nd2))], type = "l", lty = 1,
    col = redblue1[latcat], lwd = 2, axes = FALSE,
    xlab = "Time", ylab = expression(paste("Effects on ", mu)),
    ylim = ylim.mu)
  lines(fdaymu ~ day, data = nd4, col = "black", lwd = 2, lty = 2)
  legend("bottom", c("Latitude", paste(latbreaks[3:1], latbreaks[4:2], sep = "-")),
    col = c(NA, rev(redblue2)), lty = 1, lwd = 2, bty = "n", cex = 0.8)
  legend("topright", c("Mean effect", "Spatial-varying effect"),
    col = "black", lty = 2:1, box.col = NA, bg = NA, cex = 0.7) 
  abline(h = 0, lwd = 0.1)
  Axis(as.Date(c("1970-01-01", "1970-12-31")), side = 1)
  axis(2)
  box()
  main(expression(paste("Seasonal ", mu, " effect")), line = 0.3)
  dev.off()

  epng("figures/rainmodel-effects-season-sigma.png", width = 4.5, height = 3.5)
  par(mar = c(4.1, 4.1, 1.5, 1.5))
  matplot(nd2$day, nd2[, grep("fsigma", names(nd2))], type = "l", lty = 1,
    col = redblue1[latcat], lwd = 2, axes = FALSE,
    xlab = "Time", ylab = expression(paste("Effects on ", log(sigma))),
    ylim = ylim.sigma)
  lines(fdaysigma ~ day, data = nd4, col = "black", lwd = 2, lty = 2)
  abline(h = 0, lwd = 0.1)
  Axis(as.Date(c("1970-01-01", "1970-12-31")), side = 1)
  axis(2)
  box()
  main(expression(paste("Seasonal ", sigma, " effect")), line = 0.3)
  dev.off()

  epng("figures/rainmodel-effects-spatial-mu.png", width = 4.5, height = 3.5)
  par(mar = c(4.1, 4.1, 1.5, 1.5))
  plot(Austria, xlab = "Longitude [deg]", ylab = "Latitude [deg]", lwd = 0.3)
  bamlss:::xymap(lon, lat, fmu, data = nd, pos = "topleft", layout = FALSE, map = FALSE,
    add = TRUE, color = diverge_hcl, shift = c(0.09, 0.04), distance.labels = 0,
    width = 0.3, symmetric = TRUE, swap = FALSE, range = ylim.mu, digits = 1)
  plot(Austria, add = TRUE, lwd = 0.3)
  plotUnionAustria(add = TRUE)
  box()
  axis(1)
  axis(2)
  main(expression(paste("Spatial ", mu, " effect")), line = 0.3)
  dev.off()

  epng("figures/rainmodel-effects-spatial-sigma.png", width = 4.5, height = 3.5)
  par(mar = c(4.1, 4.1, 1.5, 1.5))
  plot(Austria, xlab = "Longitude [deg]", ylab = "Latitude [deg]", lwd = 0.3)
  bamlss:::xymap(lon, lat, fsigma, data = nd, pos = "topleft", layout = FALSE, map = FALSE,
    add = TRUE, color = diverge_hcl, shift = c(0.09, 0.04), distance.labels = 0,
    width = 0.3, symmetric = TRUE, swap = FALSE, range = ylim.sigma, digits = 1)
  plot(Austria, add = TRUE, lwd = 0.3)
  plotUnionAustria(add = TRUE)
  box()
  axis(1)
  axis(2)
  box()
  main(expression(paste("Spatial ", sigma, " effect")), line = 0.3)
  dev.off()

  epng("figures/rainmodel-effects-predict.png", width = 9, height = 5, res = 200)
  par(mar = c(4.1, 4.1, 1.5, 1.5))
  plot(Austria, xlab = "Longitude [deg]", ylab = "Latitude [deg]", lwd = 0.3)
  bamlss:::xymap(lon, lat, rain10, data = nd, pos = "topleft", layout = FALSE, map = FALSE,
    add = TRUE, color = colors.rain, shift = c(0.1, 0.1), distance.labels = 0,
    width = 0.3 / 2, height = 0.06 * 4/6, symmetric = FALSE, swap = FALSE, digits = 1, range = c(0, 0.6),
    lrange = round(c(0, max(nd$rain10)), 1))
  plot(Austria, add = TRUE, lwd = 0.3)
  plotUnionAustria(add = TRUE)
  box()
  axis(1)
  axis(2)
  box()
  main("Mean precipitation [mm] January 10th")
  dev.off()

  plot.mean.fit <- function(stations = 1, raw = TRUE, mean = TRUE, map = FALSE,
    col.raw = rgb(0.1, 0.1, 0.1, alpha = 0.2), col.mean = rgb(1, 0, 0, alpha = 0.2))
  {
    ylim <- NULL
    for(i in stations)
      ylim <- c(ylim , range(nd3[[i]][, 2:3]))
    ylim <- range(ylim)
    if(map)
      par(mfrow = c(2, 1))
    plot(1, 1, type = "n", ylim = ylim, xlim = c(1, 365), axes = FALSE,
      xlab = "Time", ylab = "Raw vs. fitted values")
    Axis(as.Date(c("1970-01-01", "1970-12-31")), side = 1)
    axis(2)
    box()
    if(raw) {
      for(i in stations) {
        lines(mean ~ day, data = nd3[[i]], col = col.raw, lwd = 2)
      }
    }
    if(mean) {
      for(i in stations) {
        lines(pred ~ day, data = nd3[[i]], col = col.mean, lwd = 2)
      }
    }
    if(map) {
      plot(Austria)
      for(i in stations) {
        co <- attr(nd3[[i]], "co")
        points(co$lon, co$lat, pch = 16, cex = 1.5)
        text(co$lon, co$lat, i, pos = 3, cex = 1.5)
      }
    }
  }
}


###################################
## (2) Fire response time model. ##
###################################
if(!file.exists("figures/firemodel-max-acf.png")) {
  library("bamlss")
  library("sp")
  library("maptools")
  library("rgeos")

  data("LondonFire", package = "bamlss")
  load("firemodel_plotdata.rda")

  plot.daytime <- function(x, daytime = 1, n = 30, range = NULL, lrange = NULL, main = NULL, cores = 4, chunks = 1, ...)
  {
    options(warn = -1)
    gpclibPermit()
    xy <- bbox(LondonBoroughs)
    co <- expand.grid(
      "lon" = seq(min(xy[1, 1]), max(xy[1, 2]), length = n),
      "lat" = seq(min(xy[2, 1]), max(xy[2, 2]), length = n)
    )
    ob <- unionSpatialPolygons(LondonBoroughs,
      rep(1L, length = length(LondonBoroughs)))
    nob <- length(slot(slot(ob, "polygons")[[1]], "Polygons"))
    pip <- NULL
    for(j in 1:nob) {
     oco <- slot(slot(slot(ob, "polygons")[[1]], "Polygons")[[j]], "coords")
     pip <- cbind(pip, point.in.polygon(co$lon, co$lat, oco[, 1L], oco[, 2L], mode.checked = FALSE) < 1L)
    }
    pip <- apply(pip, 1, function(x) all(x))
    co <- co[pip < 1, , drop = FALSE]
    co$daytime <- daytime
    co$spatial_daytime <- predict(x, newdata = co, model = "gamma",
      term = "daytime", intercept = FALSE, cores = cores, chunks = chunks)
    plot(LondonBoroughs, xlab = "Longitude [deg]", ylab = "Latitude [deg]",
      main = NULL)
    if(is.null(lrange)) {
      lrange <- range(co$spatial_daytime)
      lrange <- c(-1 * max(abs(lrange)), max(abs(lrange)))
    }
    if(is.null(range)) {
      range <- quantile(co$spatial_daytime, probs = 0.9)
      range <- c(-1 * range, range)
    }
    bamlss:::xymap(lon, lat, spatial_daytime, data = co, pos = "bottomright",
      layout = FALSE, map = FALSE, color = diverge_hcl, swap = TRUE,
      shift = c(0.03, 0.05), symmetric = TRUE, add = TRUE,
      side.legend = 2, digits = 1, range = range, lrange = round(lrange, 1),
      width = 0.2, height = 0.04, ...)
    plot(LondonBoroughs, add = TRUE)
    box()
    axis(1)
    axis(2)

    if(!is.null(main)) main(main)
    options(warn = 0)
    return(invisible(NULL))
  }

  plot.firemodel <- function(data, what = c("curves", "spatial_prob", "spatial_td", "spatial_tc",
    "fsintens", "daytime"), spar = FALSE, main = FALSE, ...)
  {
    if(spar)
      par(mfrow = n2mfrow(length(what)), mar = c(4.1, 4.1, 3.1, 1.1))
    if("curves" %in% what) {
      matplot(data$curves[, 1], data$curves[, -1], type = "l", lty = 1,
        col = rgb(0.1, 0.1, 0,1, alpha = 0.01), xlab = "Arrivaltime [min]",
        ylab = "Effect on log relative risk")
      plot2d(firemodel$results$lambda$s.effects[["ti(arrivaltime)"]], c.select = c(1, 3),
        col.lines = rainbow_hcl(1), add = TRUE, rug = FALSE, lwd = 2)
      abline(v = data$target, col = "blue", lty = 2)
      legend("bottomright", c("Mean baseline", "Spatial-varying baseline"), lwd = c(2, 1),
        col = c(rainbow_hcl(1), "black"), box.col = NA, bg = NA, cex = 0.95)
      if(main) main("Baseline-hazard effects")
    }
    if("daytime_curves" %in% what & !is.null(data$daytime)) {
      ylim <- range(data$daytime[, -1])
      ylim[1] <- ylim[1] - abs(diff(ylim)) * 0.1
      matplot(seq(0, 24, length = 100), data$daytime[, -1], type = "l", lty = 1,
        col = rgb(0.1, 0.1, 0,1, alpha = 0.01), xlab = "Time of day [h]",
        ylab = "Effect of time of day", ylim = ylim, axes = FALSE)
      plot2d(firemodel$results$gamma$s.effects[["ti(daytime)"]], c.select = c(1, 3),
        col.lines = rainbow_hcl(1), add = TRUE, rug = FALSE, lwd = 2, axes = FALSE)
      axis(1, at = c(0, 6, 12, 18, 24))
      axis(2)
      box()
      abline(v = 8.5, col = "blue", lty = 2)
      legend("bottomright", c("Mean daytime effect", "Spatial-varying effect"), lwd = c(2, 1),
        col = c(rainbow_hcl(1), "black"), box.col = NA, bg = NA, cex = 0.95)
      if(main) main("Spatial-varying time of day effect")
    }
    if("spatial_prob" %in% what) {
      plot(LondonBoundaries, xlab = "Longitude [deg]", ylab = "Latitude [deg]")
      lr <- c(0, 1)
      xr <- c(0, 1)
      bamlss:::xymap(lon, lat, spatial_prob, data = data$spatial, pos = "bottomright",
        layout = FALSE, map = FALSE, color = heat_hcl,
        shift = c(0.03, 0.05), symmetric = FALSE, add = TRUE,
        side.legend = 2, digits = 1, range = xr, lrange = round(lr, 1),
        width = 0.2, height = 0.04, ...)
      plot(LondonBoroughs, add = TRUE, lwd = 0.3)
      plot(LondonBoundaries, col = NA, add = TRUE)
      box()
      axis(1)
      axis(2)
      if(main) main(paste("Prob(T > ", data$target, "); 8:30 AM", sep = ""))
    }
    if("fsintens" %in% what) {
      plot(firemodel, model = "gamma", term = "(fsintens)", spar = FALSE,
        xlab = "Fire station intensity",
        ylab = "Effect on log relative risk", rug = FALSE, scheme = 2, grid = 100)
      if(main) main("Effect of fire station intensity")
    }
    if("spatial_td" %in% what) {
      plot(LondonBoundaries, xlab = "Longitude [deg]", ylab = "Latitude [deg]")
      lr <- range(data$spatial$spatial_td)
      lr <- c(-1 * max(abs(lr)), max(abs(lr)))
      xr <- quantile(data$spatial$spatial_td, probs = 0.9)
      xr <- c(-1 * xr, xr)
      bamlss:::xymap(lon, lat, spatial_td, data = data$spatial, pos = "bottomright",
        layout = FALSE, map = FALSE, color = diverge_hcl, swap = TRUE,
        shift = c(0.03, 0.05), symmetric = TRUE, add = TRUE,
        side.legend = 2, digits = 1, range = xr, lrange = round(lr, 1),
        width = 0.2, height = 0.04, ...)
      plot(LondonBoroughs, add = TRUE, lwd = 0.3)
      plot(LondonBoundaries, col = NA, add = TRUE)
      box()
      axis(1)
      axis(2)
      if(main) main(paste("Time-varying spatial effect (t = ", data$target, ")", sep = ""))
    }
    if("daytime" %in% what) {
      plot(firemodel, model = "gamma", term = "(daytime)", spar = FALSE, rug = FALSE,
        xlab = "Time of day [h]", ylab = "Effect on log relative risk",
        scheme = 2, grid = 100, axes = FALSE)
      axis(1, at = c(0, 6, 12, 18, 24))
      axis(2)
      box()
      if(main) main("Effect of time of day")
    }
    if("spatial_tc" %in% what) {
      plot(LondonBoundaries, xlab = "Longitude [deg]", ylab = "Latitude [deg]")
      lr <- range(data$spatial$spatial_tc)
      lr <- c(-1 * max(abs(lr)), max(abs(lr)))
      xr <- quantile(data$spatial$spatial_tc, probs = 0.9)
      xr <- c(-1 * xr, xr)
      bamlss:::xymap(lon, lat, spatial_tc, data = data$spatial, pos = "bottomright",
        layout = FALSE, map = FALSE, color = diverge_hcl, swap = TRUE,
        shift = c(0.03, 0.05), symmetric = TRUE, add = TRUE,
        side.legend = 2, digits = 1, range = xr, lrange = round(lr, 1),
        width = 0.2, height = 0.04, ...)
      plot(LondonBoroughs, add = TRUE, lwd = 0.3)
      plot(LondonBoundaries, col = NA, add = TRUE)
      box()
      axis(1)
      axis(2)
      if(main) main("Time-constant spatial effect")
    }
    if("spatial_daytime" %in% what) {
      plot(LondonBoundaries, xlab = "Longitude [deg]", ylab = "Latitude [deg]")
      lr <- range(data$spatial$spatial_daytime)
      lr <- c(-1 * max(abs(lr)), max(abs(lr)))
      xr <- quantile(data$spatial$spatial_daytime, probs = 0.9)
      xr <- c(-1 * xr, xr)
      bamlss:::xymap(lon, lat, spatial_daytime, data = data$spatial, pos = "bottomright",
        layout = FALSE, map = FALSE, color = diverge_hcl, swap = TRUE,
        shift = c(0.03, 0.05), symmetric = TRUE, add = TRUE,
        side.legend = 2, digits = 1, range = xr, lrange = round(lr, 1),
        width = 0.2, height = 0.04, ...)
      plot(LondonBoroughs, add = TRUE, lwd = 0.3)
      plot(LondonBoundaries, col = NA, add = TRUE)
      box()
      axis(1)
      axis(2)
      if(main) main("Spatial daytime effect")
    }
    if("stations" %in% what) {
      co <- bbox(LondonBoroughs)
      scale <- 0.5
      ylim <- c(co["y", 1], co["y", 2] + scale * abs(diff(co["y", ])))
      plot(LondonBoroughs, main = "",
        xlab = "Longitude [deg]", ylab = "", ylim = ylim,
        col = gray(0.9), lwd = 0.3)
      cfun <- function(n) { heat_hcl(n, alpha = 0.8) }
      pal <- make_pal(col = cfun, ncol = 99, data = LondonFire$arrivaltime, symmetric = FALSE,
        range = c(3, 8), swap = TRUE)
      ld <- as.data.frame(LondonFire)
      ld <- ld[order(ld$arrivaltime, decreasing = TRUE), ]
      points(ld$lon, ld$lat, pch = 4, col = pal$map(ld$arrivaltime), cex = 0.6)
      points(LondonFStations, pch = 16, col = "blue", cex = 0.8)
      points(LondonFStations, cex = 0.8)
      plot(LondonBoundaries, col = NA, add = TRUE)
      box()
      axis(1)
      axis(4, at = round(seq(co[2, 1] + 0.07, co[2, 2] - 0.07, length = 4), 2))
      mtext(paste("Latitude [deg]",
        paste(rep(" ", 30), collapse = "", sep = "")), side = 4, line = 3)
      legend("bottomleft", c("Fire", "Station"), pch = c(4, 21),
        col = c("black", "black"), pt.bg = c(NA, "blue"), box.col = NA, bg = NA, cex = 0.95)
      i <- order(LondonFire$arrivaltime, decreasing = TRUE)
      atimes <- LondonFire$arrivaltime[i]
      lon <- LondonFire$lon[i]
      shift <- scale * 0.1 * abs(diff(co["y", ]))
      y <- scale2(atimes, co["y", 2] + shift, ylim[2])
      rect(lon - 0.0015, rep(co["y", 2] + shift,length = length(y)), lon + 0.0015, y,
        col = pal$map(atimes), border = NA)
      lines(co[1, ], rep(co["y", 2] + shift, 2))
      wt <- c(min(atimes), 6, 12, 20)
      wt2 <- scale2(wt, co["y", 2] + shift, ylim[2])
      axis(2, at = wt2, labels = fmt(wt, 2, 0))
      if(main) main(paste(paste(rep(" ", 70), collapse = "", sep = ""), "Arrivaltime [min]"))
    }
  }

  plot.griddata <- function(n = 800, FUN = NULL,
    color = heat_hcl, symmetric = FALSE, swap = TRUE, type = c("hexagonal", "regular"),
    main = NULL, xlab = "Longitude [deg]", ylab = "Latitude [deg]", variable = "arrivaltime", ...)
  {
    if(is.null(FUN))
      FUN <- function(x) { sum(x > 6, na.rm = TRUE) / length(x) }

    type <- match.arg(type)

    LBP0 <- unionSpatialPolygons(LondonBoroughs,
      rep(1L, length = length(LondonBoroughs)))
    LBP <- spTransform(LBP0, CRS("+init=epsg:27700"))
    xy <- bbox(LondonBoroughs)
    dx <- abs(diff(range(xy[1, ])))
    dy <- abs(diff(range(xy[2, ])))
    scale <- 0.05
    co <- expand.grid(
      "lon" = seq(min(xy[1, 1]) - scale * dx, max(xy[1, 2]) + scale * dx, length = n),
      "lat" = seq(min(xy[2, 1]) - scale * dy, max(xy[2, 2]) + scale * dy, length = n)
    )
    coordinates(co) <- c("lon", "lat")
    proj4string(co) <- CRS("+init=epsg:4326")
    co <- spTransform(co, CRS("+init=epsg:27700"))
    Pts <- spsample(co, n = n, type = type, offset = c(0.5, 0.5))
    if(type == "hexagonal") {
      Pols <- HexPoints2SpatialPolygons(Pts)
    } else {
      Pols <- points2grid(Pts)
      Pols <- as(Pols, "SpatialPolygons")
      proj4string(Pols) <- CRS("+init=epsg:27700")
    }
    Pols <- spTransform(Pols, CRS("+init=epsg:4326"))

    clip <- gIntersection(LBP0, Pols, byid = TRUE, drop_lower_td = TRUE)
    agg <- aggregate(LondonFire[variable], clip, FUN = FUN)
    
    colors <- colorlegend(x = agg[[variable]], plot = FALSE,
      color = color, symmetric = symmetric, swap = swap, ...)

    plot(LondonBoundaries, col = gray(0.9), xlab = xlab, ylab = ylab,
      main = main)
    plot(agg, col = colors$map(agg[[variable]]), add = TRUE, border = NA)
    plot(LondonBoundaries, col = NA, add = TRUE)
    plot(LondonBoroughs, lwd = 0.3, add = TRUE)

    colorlegend(x = agg[[variable]], plot = FALSE, add = TRUE,
      color = color, swap = swap, symmetric = symmetric, width = 0.2, height = 0.04,
      pos = "bottomright", side.legend = 2, shift = c(0.03, 0.05), ...)
  }

  plot.firedata <- function()
  {
    par(mfrow = c(2, 2), mar = rep(0, 4), oma = c(4.1, 4.1, 4.1, 4.1))
    plot(LondonBoundaries, col = gray(0.9))
    plot(LondonBoroughs, lwd = 0.3, add = TRUE)
    points(LondonFire, pch = 4, col = rgb(1, 0, 0, alpha = 0.3), cex = 0.4)
    points(LondonFStations, pch = 16, col = rgb(0, 0, 1), cex = 0.8)
    points(LondonFStations, cex = 0.8)
    box()
    axis(2)
    axis(3)
    legend("topright", c("Fire", "Station"), pch = c(4, 21),
      col = c(rgb(1, 0, 0, alpha = 1), "black"),
      pt.bg = c(NA, rgb(0, 0, 1)), box.col = NA, bg = NA, cex = 0.95)

    plot.griddata(type = "hexagonal", at = c(0, 0.5, 1), labels = c("0.0", "0.5", "1.0"),
      main = "", xlab = "", ylab = "")
    box()
    axis(3)
    legend("topright", "% Arrivatime > 6min", box.col = NA, bg = NA, cex = 0.95)

    plot.griddata(type = "hexagonal", digits = 0,
      main = "", xlab = "", ylab = "", FUN = length, range = c(1, 60), lrange = c(1, 80))
    box()
    axis(2)
    legend("topright", "#Number of fires", box.col = NA, bg = NA, cex = 0.95)

    hist(LondonFire$arrivaltime, breaks = 50, col = gray(0.9), freq = FALSE,
      xlab = "", main = "", axes = FALSE)
    abline(v = 6, lty = 2, col = "blue", lwd = 2)
    box()
    axis(1)
    axis(4)
    mtext("Arrivaltime [min]", side = 1, line = 3)
    mtext("Density", side = 4, line = 3)
    F <- ecdf(LondonFire$arrivaltime)
    P <- round(1 - F(6), 2)
    legend("topright", "t = 6min",
      lwd = 2, lty = 2, col = "blue",
      title = paste("Prob(T > 6min) =", P),
      box.col = NA, bg = NA, cex = 0.95)

    mtext("Latitude [deg]", side = 2, line = 3, outer = TRUE)
    mtext("Longitude [deg]", side = 3, line = 3, outer = TRUE)
  }

  png("figures/firemodel-data.png", units = "in", res = 120, width = 8 * 1.25, height = 6.5 * 1.25)
  plot.firedata()
  dev.off()

  epng("figures/firemodel-effects-baseline.png", width = 4.5, height = 3.5)
  par(mar = c(4.1, 4.1, 1.5, 1.5))
  plot.firemodel(firemodel_plotdata, what = "curves", main = TRUE)
  dev.off()

  epng("figures/firemodel-effects-fsintens.png", width = 4.5, height = 3.5)
  par(mar = c(4.1, 4.1, 1.5, 1.5))
  plot.firemodel(firemodel_plotdata, what = "fsintens", main = TRUE)
  dev.off()

  epng("figures/firemodel-effects-daytime.png", width = 4.5, height = 3.5)
  par(mar = c(4.1, 4.1, 1.5, 1.5))
  plot.firemodel(firemodel_plotdata, what = "daytime", main = TRUE)
  dev.off()

  epng("figures/firemodel-effects-daytime-curves.png", width = 4.5, height = 3.5)
  par(mar = c(4.1, 4.1, 1.5, 1.5))
  plot.firemodel(firemodel_plotdata, what = "daytime_curves", main = TRUE)
  dev.off()

  epng("figures/firemodel-effects-prob.png", width = 4.5, height = 3.5)
  par(mar = c(4.1, 4.1, 1.5, 1.5))
  plot.firemodel(firemodel_plotdata, what = "spatial_prob", main = TRUE)
  dev.off()

  epng("figures/firemodel-effects-spatial-td.png", width = 4.5, height = 3.5)
  par(mar = c(4.1, 4.1, 1.5, 1.5))
  plot.firemodel(firemodel_plotdata, what = "spatial_td", main = TRUE)
  dev.off()

  epng("figures/firemodel-effects-spatial-tc.png", width = 4.5, height = 3.5)
  par(mar = c(4.1, 4.1, 1.5, 1.5))
  plot.firemodel(firemodel_plotdata, what = "spatial_tc", main = TRUE)
  dev.off()

  epng("figures/firemodel-max-acf.png", width = 4.5, height = 3.5)
  par(mar = c(4.1, 4.1, 1.5, 1.5))
  plot(firemodel, which = "max-acf", thin = 4, lag = 100)
  dev.off()

  lr <- range(firemodel_plotdata$spatial$spatial_tc)
  lr <- c(-1 * max(abs(lr)), max(abs(lr)))
  rr <- quantile(firemodel_plotdata$spatial$spatial_tc, probs = 0.5)
  rr <- round(c(-1 * rr, rr), 1)

  target <- c(0, 3, 6, 9, 12, 15, 18, 21, 24)
  tmain <- c("00:00", "03:00", "06:00", "09:00", "12:00", "15:00", "18:00", "21:00", "24:00")
  tmain <- paste("Time of day", tmain)
  firemodel$family <- cox_bamlss()
  for(i in seq_along(target)) {
    cat("create figure for target", i, "\n")
    epng(paste("figures/firemodel-daytime-t", i, ".png", sep = ""), width = 4.5, height = 3.5)
    par(mar = c(4.1, 4.1, 1.5, 1.5))
    plot.daytime(firemodel, daytime = target[i], n = 120, range = rr, lrange = round(lr, 1), main = tmain[i])
    dev.off()
  }
}

