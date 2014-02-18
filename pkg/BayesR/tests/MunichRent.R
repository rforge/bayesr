library("BayesR")

data("rent99", "MunichBnd", package = "BayesR")
rent99$rent <- rent99$rent / 1000

f <- list(
  rent ~ bath + kitchen + location + heating +
    sx(area) + sx(yoc) + sx(district, bs = "mrf", map = MunichBnd),
  sigma ~ bath + kitchen + location + heating +
    sx(area) + sx(yoc) + sx(district, bs = "mrf", map = MunichBnd)
)

b1 <- bayesr(f, family = gamma, data = rent99, engine = "BayesX", verbose = TRUE)


ppos <- centroids(MunichBnd, id = rent99$district)
rent99$x <- ppos$x
rent99$y <- ppos$y

f <- list(
  rent ~ bath + kitchen + location + heating +
    s(area) + s(yoc) + s(x, y),
  sigma ~ bath + kitchen + location + heating +
    s(area) + s(yoc) + s(x, y)
)

b2 <- bayesr(f, family = gamma, data = rent99, engine = "IWLS", method = c("backfitting", "MCMC"))

nd <- centroids(MunichBnd)
nd$fmu <- predict(b2, newdata = nd, model = "mu", term = "s(x,y)", intercept = FALSE)
nd$fsigma <- predict(b2, newdata = nd, model = "sigma", term = "s(x,y)", intercept = FALSE)

nd$district <- names(MunichBnd)

par(mfrow = c(1, 2))
plotmap(MunichBnd, x = nd[, c("district", "fmu")],
  interp = TRUE, type = "mba", range = c(-0.1, 0.1),
  pos = "bottomleft", side.legend = 2, side.ticks = 2,
  width = 0.25, height = 0.05, distance.labels = 0.01)

plotmap(MunichBnd, x = nd[, c("district", "fsigma")],
  interp = TRUE, type = "mba", range = c(-0.3, 0.3),
  pos = "bottomleft", side.legend = 2, side.ticks = 2,
  width = 0.25, height = 0.05, distance.labels = 0.01)

