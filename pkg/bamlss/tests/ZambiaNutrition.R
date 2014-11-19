library("bamlss")

data("ZambiaNutrition", package = "R2BayesX")
data("ZambiaBnd", package = "R2BayesX")

K <- neighbormatrix(ZambiaBnd, id = ZambiaNutrition$district)

ZambiaNutrition$district <- ZambiaNutrition$district2 <- as.factor(ZambiaNutrition$district)

f <- stunting ~ memployment + urban + gender + meducation + s(mbmi) +
  s(agechild, k = 40) + s(district, bs = "mrf", xt = list("penalty" = K)) +
  s(district2, bs = "re")

b0 <- bamlss(f, family = gaussian, data = ZambiaNutrition, engine = "JAGS")
b1 <- bamlss(f, family = gaussian, data = ZambiaNutrition, method = "backfitting", update = "optim")

f <- stunting ~ -1 + memployment + urban + gender + meducation + sx(mbmi) +
  sx(agechild) + sx(district, bs = "mrf", map = ZambiaBnd) +
  sx(district, bs = "re")

b2 <- bamlss(f, district ~ 1, family = gaussian, data = ZambiaNutrition, engine = "BayesX")

