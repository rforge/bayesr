##############################
## (1) Precipitation model. ##
##############################
## Austria: http://www.statistik.at/web_de/klassifikationen/regionale_gliederungen/nuts_einheiten/index.html
## AustriaTopo: https://www.ngdc.noaa.gov/mgg/global/global.html
if(!file.exists("rainmodel.rda") & FALSE) {
  library("bamlss")

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
  homstart <- subset(homstart, year >= 1979)

  f <- list(
    "mu" = sqrt(raw) ~ s(elevation,k=4) + ti(day,bs="cc",k=10) + ti(long,lat,bs="tp",d=2,k=30) +
      ti(day,long,lat,bs=c("cc","tp"),d=c(1,2),k=c(8,30)),
    "sigma" = ~ s(elevation,k=4) + ti(day,bs="cc",k=10) + ti(long,lat,bs="tp",d=2,k=30) +
      ti(day,long,lat,bs=c("cc","tp"),d=c(1,2),k=c(8,30))
  )

  rainmodel <- bamlss(f, data = homstart, family = "cnorm",
    binning = TRUE, before = TRUE, gam.side = FALSE, samplestats = FALSE, results = FALSE,
    eps = 0.001, n.iter = 4000, burnin = 2000, thin = 10, cores = 8)

  save(rainmodel, file = "rainmodel.rda")
}


###################################
## (2) Fire response time model. ##
###################################
if(!file.exists("firemodel.rda") & FALSE) {
  library("bamlss")
  library("spatstat")
  library("sp")
  library("maptools")
  library("raster")
  library("rgeos")

  data("LondonFire")

  f <- list(
    Surv(arrivaltime) ~ ti(arrivaltime,k=20) + ti(arrivaltime,lon,lat,d=c(1,2),k=c(5,30)),
    gamma ~ s(fsintens) + ti(daytime,bs="cc",k=30) + ti(lon,lat,k=80,d=2) +
      ti(daytime,lon,lat,bs=c("cc","cr"),d=c(1,2),k=c(10,30))
  )

  firemodel <- bamlss(f, data = LondonFire, family = "cox",
    subdivisions = 25, nu = 0.01, maxit = 1000,
    n.iter = 4000, burnin = 2000, thin = 10, cores = 8)

  save(firemodel, file = "firemodel.rda")

  data_basehaz <- function(n = 30, target = 6, k = 20, ...)
  {
    gpclibPermit()

    spatial_daytime <- grepl("(daytime,lon,lat)",
      paste(all.labels.formula(terms(firemodel, model = "gamma")), collapse = "+"), fixed = TRUE)

    firemodel$family <- cox.bamlss()

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
    nd <- NULL
    atime <- with(as.data.frame(LondonFire), seq(min(arrivaltime), max(arrivaltime), length = 100))
    if(spatial_daytime)
      dtime <- with(as.data.frame(LondonFire), seq(min(daytime), max(daytime), length = 100))
    for(i in 1:nrow(co)) {
      td <- if(spatial_daytime) {
        data.frame("arrivaltime" = atime, "daytime" = dtime)
      } else data.frame("arrivaltime" = atime)
      td$lon <- co$lon[i]
      td$lat <- co$lat[i]
      td$id <- i
      nd <- rbind(nd, td)
    }

    nd$id <- as.factor(nd$id)
    nd$p50atime <- predict(firemodel, newdata = nd, model = "lambda", FUN = mean, intercept = FALSE)
    if(spatial_daytime)
      nd$p50dtime <- predict(firemodel, newdata = nd, model = "gamma", term = "daytime", intercept = FALSE)
    nd$atime <- nd$arrivaltime
    if(spatial_daytime)
      nd$dtime <- nd$daytime

    fbh <- fdt <- NULL
    for(i in unique(nd$id)) {
      j <- nd$id == i
      fbh <- cbind(fbh, nd$p50atime[j])
      if(spatial_daytime)
        fdt <- cbind(fdt, nd$p50dtime[j])
    }
    fbh <- cbind(atime, fbh)
    if(spatial_daytime)
      fdt <- cbind(dtime, fdt)

    nd$daytime <- 8.5
    nd$arrivaltime <- target
    nd <- unique(nd[, c("lon", "lat", "arrivaltime", "daytime")])

    co <- bbox(LondonBoroughs)
    co <- owin(co["x", ], co["y", ])
    fsppp <- ppp(x = LondonFStations$lon, LondonFStations$lat, window = co)
    fsintens <- density.ppp(fsppp, bw.diggle)
    fsintens <- raster(fsintens)
    proj4string(fsintens) <- CRS("+init=epsg:4326")
    nd$fsintens <- extract(fsintens, as.matrix(nd[ , c("lon", "lat")]))

    nd$spatial_prob <- predict(firemodel, newdata = nd,
      term = c("(arrivaltime)", "(arrivaltime,lon,lat)", "(fsintens)", "(daytime)", "(lon,lat)"),
      intercept = TRUE, type = "prob", time = target, ...)
    nd$spatial_tc <- predict(firemodel, newdata = nd, model = "gamma",
      term = "(lon,lat)", intercept = FALSE)
    nd$spatial_td <- predict(firemodel, newdata = nd, model = "lambda",
      term = "(arrivaltime,lon,lat)", intercept = FALSE)
    if(spatial_daytime) {
      nd$spatial_daytime <- predict(firemodel, newdata = nd, model = "gamma",
        term = "(daytime,lon,lat)", intercept = FALSE)
    }

    return(list("curves" = fbh, "daytime" = fdt, "spatial" = nd, "target" = target))
  }

  firemodel_plotdata <- data_basehaz(120, 6, subdivisions = 25)

  save(firemodel, firemodel_plotdata, file = "firemodel_plotdata.rda")
}

