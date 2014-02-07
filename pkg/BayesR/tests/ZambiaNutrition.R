library("BayesR")

data("ZambiaNutrition", package = "R2BayesX")
data("ZambiaBnd", package = "R2BayesX")

K <- neighbormatrix(ZambiaBnd, id = ZambiaNutrition$district)

ZambiaNutrition$district <- ZambiaNutrition$district2 <- as.factor(ZambiaNutrition$district)

f <- stunting ~ memployment + urban + gender + meducation + s(mbmi) +
  s(agechild, k = 40) + s(district, bs = "mrf", xt = list("penalty" = K)) +
  s(district2, bs = "re")

b0 <- bayesr(f, family = gaussian, data = ZambiaNutrition, engine = "JAGS")
b1 <- bayesr(f, f, family = gaussian, data = ZambiaNutrition,
  engine = "IWLS", method = "backfitting")

f <- stunting ~ -1 + memployment + urban + gender + meducation + sx(mbmi) +
  sx(agechild) + sx(district, bs = "mrf", map = ZambiaBnd) +
  sx(district, bs = "re")

b2 <- bayesr(f, district ~ 1, family = gaussian, data = ZambiaNutrition, engine = "BayesX")

