library("BayesR")

data("rent99", "MunichBnd", package = "BayesR")
rent99$rent <- rent99$rent / 1000

f <- rent ~ bath + kitchen + location + heating +
  sx(area) + sx(yoc) + sx(district, bs = "mrf", map = MunichBnd)

b1 <- bayesx2(f, family = gamma.BayesR, data = rent99)
b2 <- bayesx2(f, family = gaussian.BayesR, data = rent99)
