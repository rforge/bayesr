library("bamlss")

if(!file.exists("rainmodel.rda")) {
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
  homstart2 <- subset(homstart, year >= 1979)

  f <- list(
    "mu" = sqrt(raw) ~ elevation + ti(day,bs="cc") + ti(long,lat,bs="tp",d=2,k=50) +
      ti(day,long,lat,bs=c("cc","tp"),d=c(1,2),mp=FALSE,k=c(8,30)),
    "sigma" = ~ elevation + ti(day,bs="cc") + ti(long,lat,bs="tp",d=2,k=50) +
      ti(day,long,lat,bs=c("cc","tp"),d=c(1,2),mp=FALSE,k=c(8,30))
  )

  rainmodel <- bamlss(f, data = homstart2, family = "cnorm",
    binning = TRUE, before = TRUE, gam.side = FALSE,
    samplestats = FALSE, results = FALSE,
    eps = 0.001, n.iter = 6000, burnin = 4000, thin = 10, cores = 7)

  save(rainmodel, file = "rainmodel.rda")
}
