library("bamlss")

data("india", "india.bnd", package = "gamboostLSS")

K <- neighbormatrix(india.bnd, id = india$mcdist)

f <- list(
  stunting ~ s(cbmi) + s(cage) + s(mbmi) + s(mage) + s(mcdist,bs="mrf",xt=list(penalty=K)),
  sigma ~ s(cbmi) + s(cage) + s(mbmi) + s(mage) + s(mcdist,bs="mrf",xt=list(penalty=K))
)

b <- bamlss(f, data = india)

nd <- data.frame("mcdist" = unique(india$mcdist))

nd$sp.mu <- predict(b, newdata = nd, model = "mu", term = "s(mcdist)", intercept = FALSE)
nd$sp.sigma <- predict(b, newdata = nd, model = "sigma", term = "s(mcdist)", intercept = FALSE)

par(mfrow = c(1, 2), mar = rep(0, 4))
plotmap(india.bnd, x = nd$sp.mu, id = nd$mcdist,
  pos = "bottomright", color = diverge_hcl, symmetric = TRUE,
  side.legend = 2, swap = FALSE, shift = c(0.2, 0.1))
plotmap(india.bnd, x = nd$sp.sigma, id = nd$mcdist,
  pos = "bottomright", color = diverge_hcl, symmetric = TRUE,
  side.legend = 2, swap = FALSE, shift = c(0.2, 0.1))
