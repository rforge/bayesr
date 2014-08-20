library("BayesR")

data("india", "india.bnd", package = "gamboostLSS")
india <- cbind(india, centroids(india.bnd, id = india$mcdist))
india$stunting2 <- scale(india$stunting)

f <- list(
  stunting2 ~ s(cbmi) + s(cage) + s(mbmi) + s(mage) + s(x, y),
  sigma2 ~ s(cbmi) + s(cage) + s(mbmi) + s(mage) + s(x, y)
)

b <- bayesr(f, data = india, method = c("backfitting", "MCMC"),
  update = "iwls", propose = "iwls", inner = TRUE,
  n.iter = 1200, burnin = 200, thin = 1)

india$mu.sp <- predict(b, model = "mu", term = "s(x,y)", intercept = FALSE)
india$sigma2.sp <- predict(b, model = "sigma2", term = "s(x,y)", intercept = FALSE)

plotmap(india.bnd, x = india$mu.sp, id = india$mcdist)
plotmap(india.bnd, x = india$sigma2.sp, id = india$mcdist)
