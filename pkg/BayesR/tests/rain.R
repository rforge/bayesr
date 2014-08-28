library("BayesR")

dpath <- system.file(package = "BayesR", "data")
if(!file.exists(file <- file.path(dpath, "homstart.rda"))) {
  homstart_data(dir = file, load = TRUE)
} else load(file)

rain <- na.omit(subset(homstart, year >= 2000)[, c("raw", "cens", "long", "lat", "day", "elevation")])
rain$raw2 <- sqrt(rain$raw)
rain$raw2[is.na(rain$raw2)] <- 0

f <- list(
  sqrt(cens) ~ s(day, bs = "cc") + s(elevation) + s(long, lat),
  sigma2 ~ s(day, bs = "cc") + s(elevation) + s(long, lat)
)

b0 <- bayesr(f, data = rain, method = "backfitting", update = "iwls", sample = "iwls")

f <- list(
  raw ~ s(day, bs = "cc") + s(elevation),
  sigma ~ s(day, bs = "cc") + s(elevation)
)

b1 <- bayesr(f, data = rain, family = gF(tobit), method = "MP", n.samples = 10)
