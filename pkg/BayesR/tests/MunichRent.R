library("bamlss")

data("rent99", "MunichBnd", package = "bamlss")
rent99$rent <- rent99$rent / 1000

f <- list(
  rent ~ bath + kitchen + location + cheating +
    sx(area) + sx(yearc) + sx(district, bs = "mrf", map = MunichBnd),
  sigma ~ bath + kitchen + location + cheating +
    sx(area) + sx(yearc) + sx(district, bs = "mrf", map = MunichBnd)
)

b1 <- bamlss(f, family = gamma, data = rent99, engine = "BayesX", verbose = TRUE)

rent99 <- cbind(rent99, centroids(MunichBnd, id = rent99$district))

f <- list(
  rent ~ bath + kitchen + location + cheating +
    s(area) + s(yearc) + s(x, y, k = 100),
  sigma ~ bath + kitchen + location + cheating +
    s(area) + s(yearc) + s(x, y, k = 100)
)

b2 <- bamlss(f, family = gamma, data = rent99, method = c("MP", "MCMC"), propose = "iwls")


nd <- centroids(MunichBnd)
nd$fmu <- predict(b2, newdata = nd, model = "mu", term = "s(x,y)", intercept = FALSE)
nd$fsigma <- predict(b2, newdata = nd, model = "sigma", term = "s(x,y)", intercept = FALSE)

nd$district <- names(MunichBnd)

par(mfrow = c(2, 2))
plot(b1, model = "mu", term = 3, map = MunichBnd,
  range = c(-0.1, 0.1), interp = TRUE, type = "mba",
  pos = "bottomleft", side.legend = 2, side.ticks = 2,
  width = 0.25, height = 0.05, distance.labels = 0.01,
  spar = FALSE)

plot(b1, model = "sigma", term = 3, map = MunichBnd,
  range = c(-0.3, 0.3), interp = TRUE, type = "mba",
  pos = "bottomleft", side.legend = 2, side.ticks = 2,
  width = 0.25, height = 0.05, distance.labels = 0.01,
  spar = FALSE)

plotmap(MunichBnd, x = nd[, c("district", "fmu")],
  range = c(-0.1, 0.1), interp = TRUE, type = "mba",
  pos = "bottomleft", side.legend = 2, side.ticks = 2,
  width = 0.25, height = 0.05,
  distance.labels = 0.01)

plotmap(MunichBnd, x = nd[, c("district", "fsigma")],
  range = c(-0.3, 0.3), interp = TRUE, type = "mba",
  pos = "bottomleft", side.legend = 2, side.ticks = 2,
  width = 0.25, height = 0.05,
  distance.labels = 0.01)


grid <- 200
bbox <- bbox(bnd2sp(MunichBnd))
nd <- expand.grid("x" = seq(bbox["x", 1], bbox["x", 2], length = grid),
  "y" = seq(bbox["y", 1], bbox["y", 2], length = grid))

nd$fmu <- predict(b2, newdata = nd, model = "mu", term = "s(x,y)")

i <- drop2poly(nd$x, nd$y, MunichBnd)
nd <- nd[i, ]

par(mfrow = c(1, 2), mar = rep(0, 4))
xymap(x, y, fmu, data = nd, col = heat_hcl, symmetric = FALSE,
  swap = TRUE, pos = "bottomleft", range = c(0.5, 0.59),
  side.legend = 2, side.ticks = 2, width = 0.25, height = 0.03,
  distance.labels = 0.01, layout = FALSE)
plot(MunichBnd, add = TRUE)

xymap(x, y, fsigma, data = nd, col = diverge_hcl, symmetric = TRUE,
  range = c(-0.3, 0.3), swap = FALSE, pos = "bottomleft",
  side.legend = 2, side.ticks = 2, width = 0.25, height = 0.03,
  distance.labels = 0.01, layout = FALSE)
plot(MunichBnd, add = TRUE)

