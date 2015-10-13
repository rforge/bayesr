library("bamlss")
data("GAMart", package = "R2BayesX")

## Model formulae.
bamlss.formula(num ~ s(x1) + id)

f <- list(
  num ~ s(x1) + s(x3) + id,
  sigma ~ s(x2),
  id ~ s(x3) + long,
  long ~ te(x1, x2)
)

bamlss.formula(f)
bamlss.formula(f, family = gamma.bamlss())


## Automatic filling with intercept.
bamlss.formula(num ~ x1, family = zinb.bamlss())
bamlss.formula(f, family = zinb.bamlss())


## Terms object from bamlss.formula().
terms(bamlss.formula(f))
terms.bamlss(f)
terms(bamlss.formula(f), model = c(1, 1), sterms = FALSE)
terms(bamlss.formula(f), model = c(1, 1), sterms = FALSE, pterms = FALSE)


## Get the bamlss frame.
bf <- bamlss.frame(f, data = GAMart, family = "gaussian")
names(bf)
print(bf)
head(model.frame(bf))
bf$model.frame <- NULL
head(model.frame(bf))
model.response(model.frame(bf))
response.name(bf)
response.name(model.frame(bf))


## model.matrix() and smooth.construct()
model.matrix(bf)
head(model.matrix(bf, model = c(1, 1)))
smooth.construct(bf)
bf$model.frame <- NULL
head(model.matrix(bf, model = c(1, 1)))
str(smooth.construct(bf))
str(smooth.construct(bf, model = c(1, 1)))

bf <- bamlss.frame(f, data = GAMart, family = "gaussian",
  model.matrix = FALSE, smooth.construct = FALSE)
print(bf)
model.matrix(bf)
head(model.matrix(bf, model = c(1, 1)))
smooth.construct(bf)
names(smooth.construct(bf, model = c(1, 1)))


## Complex multilevel structures.
f <- list(
  cat ~ s(x1) + s(x2) + id,
  id ~ s(x3)
)

bf <- bamlss.frame(f, data = GAMart, family = "multinomial")
print(bf)
head(model.response(model.frame(bf)))

terms(bf)
terms(formula(bf))
terms(bf, model = c(2, 1))
terms(bf, model = c(2, 1), sterms = FALSE)
terms(bf, model = c(2, 1), pterms = FALSE)

d <- GAMart[1:30, ]

model.matrix(formula(bf), data = d)
head(model.matrix(formula(bf), data = d, model = c(1, 1)))
head(model.matrix(formula(bf, model = c(1, 1)), data = d))
head(model.matrix(terms(bf, model = c(1, 1), drop = FALSE), data = d))

smooth.construct(formula(bf), data = d)
smooth.construct(terms(bf), data = d)
smooth.construct(formula(bf), data = d, model = c(1, 1))
smooth.construct(formula(bf, model = c(1, 1)), data = d)
smooth.construct(terms(bf, model = c(1, 1), drop = FALSE), data = d)


## Extract or initiallize parameters.
bf <- bamlss.frame(num|sigma ~ s(x1) + s(x2)|s(x3), data = GAMart)
p <- parameters(bf)
unlist(p)
unlist(parameters(randomize(bf)))


## Estimate model.
data("GAMart", package = "R2BayesX")
b <- bamlss(num|sigma ~ s(x1) + s(x2) + x3 + id | s(x1) + x2, data = GAMart, cores = 3, chains = 2)
samps <- samples(b)
samps <- samples(b, model = 1, term = 1)
head(samps)
plot(b, model = 1, term = 1, which = "samples")


## Predict.
p <- predict(b, model = "mu", term = "s(x2)")
plot2d(p ~ x2, data = GAMart)
p <- predict(b, model = "mu", term = "s(x2)", FUN = quantile, probs = c(0.025, 0.5, 0.975))
plot2d(p ~ x2, data = GAMart)
p <- predict(b, model = 1, term = "(Intercept)")
p <- predict(b)
names(p)
p <- predict(b, model = "mu", term = "id", intercept = FALSE)
plotblock(p ~ id, data = GAMart)

nd <- data.frame("x2" = seq(0, 1, length = 100))
nd$p <- predict(b, newdata = nd, model = "mu", term = "s(x2)")
plot(nd, type = "l")

nd <- data.frame("x2" = seq(0, 1, length = 100), "x3" = seq(0, 1, length = 100))
p <- predict(b, newdata = nd, model = 1, term = c("x2", "x3"), FUN = quantile)
print(head(p))

str(predict(b, term = 1:2))


## Multinomial example.
data("marital.nz", package = "VGAM")

b <- bamlss(mstatus ~ s(age), data = marital.nz,
  family = "multinomial", reference = "Married/Partnered", cores = 4)

plot(b)
plot(b, which = "samples")
plot(b, model = 1)
plot(b, model = c(1, 3))

a <- results.bamlss.default(b)
plot(a)

nd <- data.frame("age" = seq(0, 100, length = 100))
pi <- predict(b, newdata = nd, type = "parameter")
pi <- cbind(as.data.frame(pi), "MarriedPartnered" = 1)
pi <- pi / rowSums(pi)
with(nd, matplot(age, pi, type = "l", lty = 1))
legend("topright", names(pi), lwd = 1, col = 1:ncol(pi))


## Survival example.
set.seed(111)

## Sample covariates first.
n <- 600
X <- matrix(NA, nrow = n, ncol = 3)
X[, 1] <- runif(n, -1, 1)
X[, 2] <- runif(n, -3, 3)
X[, 3] <- runif(n, -1, 1)

## Specify censoring function.
cens_fct <- function(time, mean_cens) {
  ## Censoring times are independent exponentially distributed.
  censor_time <- rexp(n = length(time), rate = 1 / mean_cens)
  event <- (time <= censor_time)
  t_obs <- apply(cbind(time, censor_time), 1, min)
  ## Return matrix of observed survival times and event indicator.
  return(cbind(t_obs, event))
}

## log(time) is the baseline hazard.
lambda <-  function(time, x) {
  exp(log(time) + 0.7 * x[1] + sin(x[2]) + sin(time * 2) * x[3])
}

## Simulate data with lambda() and cens_fct().
d <- rSurvTime2(lambda, X, cens_fct, mean_cens = 5)

f <- list(
  Surv(time, event) ~ s(time, bs = "ps", k = 20) + s(time, by = x3),
  gamma ~ s(x1) + s(x2)
)

## Cox model with continuous time.
b <- bamlss(f, family = "cox", data = d, nu = 0.1, subdivisions = 50, cores = 4,
  n.iter = 4000, burnin = 1000, thin = 5)

## Predict P(T > t).
b$family <- cox.bamlss()
predict(b, type = "probabilties", time = 4)


## JAGS.
data("GAMart", package = "R2BayesX")
b <- bamlss(num|sigma ~ s(x1) + s(x2) + x3 + id | s(x1) + x2, data = GAMart, sampler = JAGS)

data("marital.nz", package = "VGAM")
b <- bamlss(mstatus ~ s(age), data = marital.nz,
  family = "multinomial", reference = "Married/Partnered",
  sampler = JAGS, cores = 4)


## Boosting.
data("GAMart", package = "R2BayesX")

f <- list(
  num ~ x1 + x2 + x3 + s(x1) + s(x2) + s(x3) + s(long, lat) + id,
  sigma ~ x1 + x2 + x3 + s(x1) + s(x2) + s(x3) + s(long, lat) + id
)

b <- bamlss(f, data = GAMart, optimizer = boost99, sampler = FALSE)
plot(b)
plot(b, which = "boost.summary")
x11()

data("india", "india.bnd", package = "gamboostLSS")

india$stunting_rs <- india$stunting / 600
K <- neighbormatrix(india.bnd, id = india$mcdist)

xt <- list("penalty" = K, "center" = FALSE, "force.center" = TRUE)

f <- list(
  stunting_rs ~ s(mage) + s(mbmi) + s(cage) + s(cbmi) + s(mcdist,bs="mrf",xt=xt),
  sigma ~ s(mage) + s(mbmi) + s(cage) + s(cbmi) + s(mcdist,bs="mrf",xt=xt)
)

b <- bamlss(f, data = india, sampler = FALSE, optimizer = boost99)


## Another survival example.
## Paper: http://arxiv.org/pdf/1503.07709.pdf
## Data: http://data.london.gov.uk/dataset/london-fire-brigade-incident-records
data("LondonFire")

f <- list(
  Surv(arrivaltime) ~ s(arrivaltime),
  gamma ~ s(fsintens, k = 20) + s(daytime, bs = "cc", k = 20) + s(lon, lat, k = 100)
)

b <- bamlss(f, data = LondonFire, family = "cox", subdivisions = 50,
  n.iter = 1200, burnin = 200, thin = 4, cores = 4)

plot(b, model = 2, term = 2, image = TRUE, grid = 200, swap = TRUE)


## Sparse matrices.
sparse.matrix <- function(n = 10, m = 20, sparse = TRUE)
{
  require("Matrix")
  x <- smooth.construct(s(x), list("x" = runif(m)), NULL)
  if(!sparse)
    return(list("X" = x$X, "S" = x$S[[1]]))
  X <- as.matrix(bdiag(rep(list(x$X), length = n)))
  S <- as.matrix(bdiag(rep(list(x$S[[1]]), length = n)))
  return(list("X" = X, "S" = S))
}

M <- sparse.matrix(sparse = TRUE)
xx <- crossprod(M$X) + 0 * M$S
b <- seq(0, 1, length = nrow(xx))
sp <- sparse.setup(M$X)
L <- sparse.chol(xx, index = list("matrix" = sp$crossprod, "ordering" = sp$ordering))
y <- sparse.forwardsolve(L, b, index = list("matrix" = sp$forward, "ordering" = sp$ordering))
z1 <- sparse.backsolve(L, y, list("matrix" = sp$backward, "ordering" = sp$ordering))
z2 <- solve(xx) %*% b
all.equal(z1, z2)
all.equal(sparse.solve(xx, b, sp), solve(xx,b))

sparse.solve(xx, diag(nrow(xx)), sp)


add2spam <- function(x, y)
{
  pointers <- attributes(x)
  x@entries <- x@entries + y[x@rowpointers, x@colindices]
  return(x)
}

