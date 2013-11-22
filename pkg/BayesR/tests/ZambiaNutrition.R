library("BayesR")

data("ZambiaNutrition", package = "R2BayesX")
data("ZambiaBnd", package = "R2BayesX")

K <- neighbormatrix(ZambiaBnd, id = ZambiaNutrition$district)

ZambiaNutrition$district <- as.factor(ZambiaNutrition$district)

f <- stunting ~ memployment + urban + gender + meducation + s(mbmi) +
  s(agechild) + s(district, bs = "mrf", xt = list("penalty" = K)) +
  s(district, bs = "re")

b0 <- bayesr(f, family = gaussian.BayesR, data = ZambiaNutrition)

f <- stunting ~ -1 + memployment + urban + gender + meducation + sx(mbmi) +
  sx(agechild) + sx(district, bs = "mrf", map = ZambiaBnd) +
  sx(district, bs = "re")

b1 <- bayesx2(f, district ~ 1, family = gaussian.BayesR, data = ZambiaNutrition)

