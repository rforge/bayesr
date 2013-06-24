dir <- path.expand("~/svn/bayesr/pkg/R2BayesX/R")
invisible(sapply(file.path(dir, list.files(dir)), source))

colorlegend(pos = "center", length.ticks = 0.1)








## multinomial example
library("R2BayesX")
data("marital.nz", package = "VGAM")
marital.nz$somefactor <- factor(sample(1:3, nrow(marital.nz), replace = TRUE),
  levels = 1:3, labels = c("red", "blue", "green"))

b <- bayesx(mstatus ~ sx(age) + somefactor, method = "MCMC",
  family = "multinomial", data = marital.nz, reference = "Married/Partnered")

fb <- predict(b, type = "response")
fb <- cbind(fb, 1 - rowSums(fb))
mycol <- c("red", "darkgreen", "blue")
mf <- model.frame(b)
ooo <- with(mf, order(age))
with(mf, matplot(age[ooo], fb[ooo, ],
  type = "l", las = 1, lwd = 2, ylim = 0:1, ylab = "Fitted probabilities",
  xlab = "Age", col = c(mycol[1], "black", mycol[-1])))
legend(x = 52.5, y = 0.62, col = c(mycol[1], "black", mycol[-1]),
  lty = 1:4, legend = levels(marital.nz$mstatus), lwd = 2)


fit.ms <- vgam(mstatus ~ s(age, df = 10), multinomial(refLevel = 2), data = marital.nz)
par(mfrow = c(3, 1))
plotvgam(fit.ms, se = TRUE, which.term = 1)
plotvgam(fit.ms, se = TRUE, which.term = 2)
plotvgam(fit.ms, se = TRUE, which.term = 3)
mu <- fitted(fit.ms)


## stepwise example
## an example of automatic model selection via null space penalization
set.seed(3); n <- 500
dat <- gamSim(1, n = n, scale = 1.2) ## simulate data
dat$x4 <- runif(n, 0, 1); dat$x5 <- runif(n, 0, 1) ## spurious

b <- gam(y ~ s(x0, bs = "ps") + s(x1, bs = "ps") + s(x2, bs = "ps") + 
  s(x3, bs = "ps") + s(x4, bs = "ps") + s(x5, bs = "ps"),
  data = dat, select = TRUE, method = "REML")

summary(b)
plot(b, pages = 1)

b2 <- bayesx(y ~ s(x0, bs = "ps") + s(x1, bs = "ps") + s(x2, bs = "ps") + 

  s(x3, bs = "ps") + s(x4, bs = "ps") + s(x5, bs = "ps"),
  data = dat, family = "gaussian", method = "STEP", dir.rm = FALSE)

b3 <- bayesx(y ~ s(x0, bs = "ps", k = 20) + s(x1, bs = "ps", k = 20) + s(x2, bs = "ps", k = 20),
  data = dat, family = "gaussian", method = "REML")


## a dummy example
## more examples
set.seed(111)
n <- 500
     
## regressors
dat <- data.frame(x = runif(n, -3, 3), z = runif(n, -3, 3),
  w = runif(n, 0, 6), fac = factor(rep(1:2, n/2)))
     
## response
dat$y <- with(dat, 1.5 + sin(x) + cos(z) * sin(w) +
  c(2.67, 3.33)[fac] + rnorm(n, sd = 0.6))

b1 <- bayesx(y ~ s(x, bs = "ps") + s(z, w, bs = "te") + fac,
  data = dat, method = "MCMC")

b2 <- bayesx(y ~ s(x, bs = "ps") + s(z, w, bs = "te") + fac,
  data = dat, method = "MCMC", contrasts = list(fac = contr.sum))


