library("BayesR")

data("rent99", "MunichBnd", package = "BayesR")
rent99$rent <- rent99$rent / 1000

f <- rent ~ bath + kitchen + location + heating +
  sx(area) + sx(yoc) + sx(district, bs = "mrf", map = MunichBnd)

b1 <- bayesr(f, family = gamma, data = rent99, engine = "BayesX")
b2 <- bayesr(f, f, family = gaussian data = rent99, engine = "BayesX")

map <- MunichBnd[unique(as.character(rent99$district))]
class(map) <- "bnd"
pen <- neighbormatrix(map)
f <- rent ~ bath + kitchen + location + heating +
  s(area) + s(yoc) + s(district, bs = "mrf", xt = list("penalty" = pen))

b3 <- bayesr(f, data = rent99, engine = "IWLS")
