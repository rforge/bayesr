library("R2BayesX")

## load zambia data and map
data("ZambiaNutrition")
data("ZambiaBnd")

## estimate model
b1 <- bayesx(stunting ~ memployment + education + urban + gender + 
  s(bmi, bs = "ps") + s(agechild, bs = "ps") +
  s(district, bs = "mrf", xt = list(map = ZambiaBnd)) + r(district),
  iter = 12000, burnin = 2000, step = 10, data = ZambiaNutrition)

## summary statistics
summary(b1)

## visualising estimation results
plot(b1, term = "s(bmi)")
plot(b1, term = "s(agechild)")

plot(b1, term = "s(district)")
plot(b1, term = "r(district)")

plot(b1, term = c("s(district)", "r(district)"), 
  map = ZambiaBnd, pos = "topleft")

## customizing graphics
plot(b1, term = "s(bmi)", resid = TRUE, cex = 0.1)
plot(b1, term = "s(agechild)", resid = TRUE, cex = 0.1)

plot(b1, term = "s(bmi)", main = "Mother body mass index", 
  xlab = "",ylab = "")
plot(b1, term = "s(bmi)", main = "Mother body mass index", 
  xlab = "",ylab = "", ylim = c(-0.8, 0.6))
plot(b1, term = "s(bmi)", main = "Mother body mass index", 
  xlab = "",ylab = "", ylim = c(-0.8, 0.6), rug = FALSE)
plot(b1, term = "s(bmi)", main = "Mother body mass index", 
  xlab = "",ylab = "", ylim = c(-0.8, 0.6), rug = FALSE, 
  col.poly = NA, lwd = 1, lty = 1)
plot(b1, term = "s(bmi)", main = "Mother body mass index", 
  xlab = "", ylab = "", ylim = c(-0.8, 0.6), rug = FALSE, 
  col.poly = NA, lwd = 1, lty = c(3, 1, 2, 2, 1))

plot(b1, term = "s(district)", map = ZambiaBnd, pos = "topleft")
plot(b1, term = "r(district)", map = ZambiaBnd, pos = "topleft")
plot(b1, term = "r(district)", map = ZambiaBnd, pos = "topleft", 
  density = 30, swap = TRUE)
plot(b1, term = "s(district)", map = ZambiaBnd, pos = "topleft", 
  names = TRUE)
plot(b1, term = "s(district)", map = ZambiaBnd, pos = "topleft", 
  names = TRUE, cex.names = 0.8, cex.legend = 0.8)
plot(b1, term = "s(district)", map = ZambiaBnd, pos = "topleft", 
  range = c(-0.2, 0.2))
plot(b1, term = "s(district)", map = ZambiaBnd, pos = "topleft", 
  c.select = "pcat95", at = c(-1, 0, 1), ncol = 3)
plot(b1, term = "s(district)", map = ZambiaBnd, pos = "topleft", 
  c.select = "pcat80", at = c(-1, 0, 1), ncol = 3)

op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
plot(b1, term = "s(district)", map = ZambiaBnd, pos = "topleft", 
  main = "Structured spatial effect", names = TRUE)
plot(b1, term = "r(district)", map = ZambiaBnd, pos = "topleft", 
  density = 40, main = "Unstructured spatial effect", names = TRUE)
par(op)

## diagnostic plots
plot(b1, which = 5:8, cex = 0.1)
plot(b1, which = "intcpt-samples")

plot(b1, term = c("memployment", "education"), 
  which = "coef-samples", ask = TRUE)

plot(b1, term = "s(bmi)", which = "coef-samples")
plot(b1, term = "s(bmi)", which = "coef-samples", acf = TRUE)

plot(b1, term = "s(bmi)", which = "var-samples")
plot(b1, term = "s(bmi)", which = "var-samples", acf = TRUE)
plot(b1, term = "s(bmi)", which = "var-samples", acf = TRUE, lag.max = 400)

plot(b1, term = "s(agechild)", which = "var-samples")
plot(b1, term = "s(district)", which = "var-samples")
plot(b1, term = "r(district)", which = "var-samples")

## sensitivity analysis 
b2 <- bayesx(stunting ~ memployment + education + urban + gender + 
  s(bmi, bs = "ps", xt = list(a = 0.00001, b = 0.00001)) +
  s(agechild, bs = "ps", xt = list(a = 0.00001, b = 0.00001)) +
  s(district, bs = "mrf", xt = list(map = ZambiaBnd, a = 0.00001, b = 0.00001)) +
  r(district, xt = list(a = 0.00001, b = 0.00001)),
  iter = 12000, burnin = 2000, step = 10, data = ZambiaNutrition)

b3 <- bayesx(stunting ~ memployment + education + urban + gender +
  s(bmi, bs = "ps", xt = list(a = 0.005, b = 0.005)) +
  s(agechild, bs = "ps", xt = list(a = 0.005, b = 0.005)) +
  s(district, bs = "mrf", xt = list(map = ZambiaBnd, a = 0.005, b = 0.005)) +
  r(district, xt = list(a = 0.005, b = 0.005)),
  iter = 12000, burnin = 2000, step = 10, data = ZambiaNutrition)

b4 <- bayesx(stunting ~ memployment + education + urban + gender +
  s(bmi, bs= "ps", xt = list(a = 0.00005, b = 0.00005)) +
  s(agechild, bs = "ps", xt = list(a = 0.00005, b = 0.00005)) +
  s(district, bs = "mrf", xt = list(map = ZambiaBnd, a = 0.00005, b = 0.00005)) +
  r(district, xt = list(a = 0.00005, b = 0.00005)),
  iter = 12000, burnin = 2000, step = 10, data = ZambiaNutrition)

plot(c(b1, b2, b3, b4), term = "s(bmi)")
plot(c(b1, b2, b3, b4), term = "s(agechild)")


## same model with REML
b5 <- bayesx(stunting ~ memployment + education + urban + gender + 
  s(bmi, bs = "ps") + s(agechild, bs = "ps") +
  s(district, bs = "mrf", xt = list(map = ZambiaBnd)) + r(district),
  method = "REML", data = ZambiaNutrition)
summary(b5)
plot(b5, map = ZambiaBnd)


## now use the stepwise algorithm
b6 <- bayesx(stunting ~ memployment + education + urban + gender + 
  s(bmi, bs = "ps") + s(agechild, bs = "ps") +
  s(district, bs = "mrf", xt = list(map = ZambiaBnd)) + r(district),
  iter = NULL, burnin = NULL, step = NULL, method = "STEP", 
  data = ZambiaNutrition, dir.rm = FALSE, CI = NULL, bootstrapsamples = NULL)
summary(b6)
