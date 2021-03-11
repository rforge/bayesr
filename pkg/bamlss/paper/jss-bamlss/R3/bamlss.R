# ---- Packages ------------------------------------------------------------
library("bamlss")
library("gamlss.dist")
library("gamlss.tr")
library("ggplot2")
library("sf")
library("rnaturalearth")
library("AER")
## install.packages("FlashAustria", repos = "http://R-Forge.R-project.org")
library("FlashAustria")


# ---- Section 2.1: Basic Bayesian regression: Logit model -----------------
data("SwissLabor", package = "AER")

## Model formula.
f <- participation ~ income + age + education + youngkids + oldkids + foreign + I(age^2)

## Estimate model.
set.seed(123)
b <- bamlss(f, family = "binomial", data = SwissLabor)

## Model summary.
summary(b)

## Compute DIC.
DIC(b)

## Figure 1: MCMC diagnostics.
plot(b, which = c("samples", "max-acf"))

## Predictions on probability scale.
nd <- data.frame(income = 11, age = seq(2, 6.2, length = 100),
  education = 12, youngkids = 1, oldkids = 1, foreign = "no")
nd$p_swiss <- predict(b, newdata = nd, type = "parameter", FUN = c95)
nd$foreign <- "yes"
nd$p_foreign <- predict(b, newdata = nd, type = "parameter", FUN = c95)

## Plot effect of age on probability.
blues <- function(n, ...) hcl.colors(n, "Blues", rev = TRUE)
plot2d(p_swiss ~ age, data = nd, ylab = "participation",
  ylim = range(c(nd$p_swiss, nd$p_foreign)), fill.select = c(0, 1, 0, 1))
plot2d(p_foreign ~ age, data = nd, add = TRUE,
  fill.select = c(0, 1, 0, 1), axes = FALSE,
  s2.col = blues, col.lines = blues(1))
legend("topright", c("Foreign", "Swiss"), lwd = 1,
  col = c(blues(1), "black"), bty = "n")


# ---- Section 2.2: Flexible model terms and estimators --------------------
## Model formula including a smooth model term for age.
f <- participation ~ income + education +
  youngkids + oldkids + foreign + s(age, k = 10)

## Estimate model.
set.seed(123)
b <- bamlss(f, family = "binomial", data = SwissLabor)

## Plot estimated smooth effect.
plot(b, term = "s(age)",
  ylab = expression(paste("Effect on Logi", t^-1, (pi))))

## Categorize age variable using quantiles.
SwissLabor$cage <- cut(SwissLabor$age,
  breaks = quantile(SwissLabor$age, prob = seq(0, 1, length = 10)),
  include.lowest = TRUE, ordered_result = TRUE)

## Model formula including the fused lasso model term for cage.
f <- participation ~ income + education + youngkids + oldkids + foreign + la(cage, fuse = 2)

## Estimate model using the lasso optimizer function.
b <- bamlss(f, family = "binomial", data = SwissLabor,
  optimizer = opt_lasso, sampler = FALSE, upper = exp(5), lower = 1,
  criterion = "AIC")

## Figure 3: BIC path, paths of coefficients and estimated effect.
par(mfrow = c(1, 3), mar = c(4.1, 4.1, 5, 7.3))
pathplot(b, 1, spar = FALSE)
pathplot(b, 2, spar = FALSE, name = "pi.s.la(cage).cage",
  main = "Coefficient paths age", lwd = 2)
page <- predict(b, term = "cage", intercept = FALSE, mstop = lasso_stop(b))
plot2d(page ~ age, data = SwissLabor, rug = TRUE, lwd = 2,
  ylab = expression(paste("Effect on Logi", t^-1, (pi))))
mtext("Estimated\nnonlinear effect", side = 3, line = 1.5, cex = 1.2, font = 2)


# ---- Section 2.3: Distributional model -----------------------------------
data("fatalities", package = "bamlss")

## Create subsets.
d19 <- subset(fatalities, year <= 2019)
d20 <- subset(fatalities, year >= 2020)

## Model formula.
f <- list(
  log(num) ~ s(week, bs = "cc", k = 20),
  sigma    ~ s(week, bs = "cc", k = 20)
)

## Estimate location-scale model.
set.seed(456)
b <- bamlss(f, data = d19, family = NO)

## Plot estimated effects.
plot(b, model = c("mu", "sigma"))

## Show predicted quantiles with 2020 data.
nd <- data.frame(week = 1:53)
par <- predict(b, newdata = nd, type = "parameter")
nd$fit <- sapply(c(0.05, 0.5, 0.95),
  FUN = function(p) { exp(family(b)$q(p, par)) })
d19w <- reshape(d19, idvar = "week",
  timevar = "year", direction = "wide")

par(mar = c(4.1, 4.1, 0.1, 0.1))
matplot(d19w$week, d19w[, -1],
  type = "o", lty = 1, pch = 16, col = gray(0.1, alpha = 0.1),
  xlab = "Week", ylab = "Number of fatalities")
grid()
matplot(nd$week, nd$fit, type = "l", lty = c(2, 1, 2),
  col = 1, lwd = 2, add = TRUE)
lines(num ~ week, data = d20, col = 2, lwd = 2, type = "o", pch = 16)



## Plot histogram of quantile residuals.
plot(b, which = "hist-resid", spar = FALSE, col = "lightgray")

## Show quantile residuals QQ-plot.
plot(b, which = "qq-resid", spar = FALSE, c95 = TRUE)


# ---- Section 4: The bamlss package ---------------------------------------
## Show optimizer/sampler information for families.
print(engines(NO, gaussian_bamlss, cox_bamlss))

## Model formula from Section 2.1.
f <- participation ~ income + age + education + youngkids + oldkids + foreign + I(age^2)

## Setup the bamlss.frame.
bf <- bamlss.frame(f, data = SwissLabor, family = "binomial")

## Run backfitting model fitting engine.
pm <- with(bf, bfit(x, y, family))

## Run MCMC engine.
set.seed(123)
samps <- with(bf, GMCMC(x, y, family, start = pm$parameters))

## Compute sampling statistics.
stats <- with(bf, samplestats(samps, x, y, family))
print(unlist(stats))


# ---- Section 4.1: The bamlss model frame ---------------------------------
## Generate simulated data set.
set.seed(111)
d <- GAMart()

## A model formula.
f <- list(
  num ~ x1 + s(x2) + s(x3) + te(lon,lat),
  sigma ~ x1 + s(x2) + s(x3) + te(lon,lat)
)

## Create the bamlss frame.
bf <- bamlss.frame(f, data = d, family = "gaussian")

## Show the structure of the bamlss frame object.
print(bf)

## Show names of smooth term list.
## Each element holds design and penalty matrices
## that can be used for estimation.
print(names(bf$x$mu$smooth.construct))


# ---- Section 4.3: Estimation engines ------------------------------------
## Show parameter names of bamlss family.
gaussian_bamlss()$names

## Example of parameter names.
paste0("mu.s.s(x3)", ".b", 1:10)


# ---- Section 5: Flexible count regression for lightning reanalysis ------
data("FlashAustria", package = "FlashAustria")

## Load estimated models.
data("FlashAustriaModel", package = "FlashAustria")

## Show first rows of the data.
head(FlashAustriaTrain)[, c("counts", "q_prof_PC1", "sqrt_cape",
  "sqrt_lsp", "d2m", "cswc_prof_PC4")]

## Number of observations.
nrow(FlashAustriaTrain)

## Generate truncate Sichel distribution.
ztSICHEL <- trun(0, family = "SICHEL", local = FALSE)

## The model formula.
f <- list(
  counts ~ s(d2m, bs = "ps") + s(q_prof_PC1, bs = "ps") +
           s(cswc_prof_PC4, bs = "ps") + s(t_prof_PC1, bs = "ps") +
           s(v_prof_PC2, bs = "ps") + s(sqrt_cape, bs = "ps"), 
         ~ s(sqrt_lsp, bs = "ps"),
         ~ 1
)

## Models are estimated with:
if(FALSE) {
  set.seed(123)
  flash_model_ztnbinom <- bamlss(f, data = FlashAustriaTrain, ## Standard interface.
    family = "ztnbinom", binning = TRUE,                      ## General arguments.
    optimizer = opt_boost, maxit = 1000,                      ## Boosting arguments.
    thin = 3, burnin = 1000, n.iter = 2000,                   ## Sampler arguments.
    light = TRUE, cores = 3)

  set.seed(123)
  flash_model_ztSICHEL <- bamlss(f, data = FlashAustriaTrain,
    family = ztSICHEL, binning = TRUE,
    optimizer = opt_boost, maxit = 1000,
    thin = 3, burnin = 1000, n.iter = 2000,
    light = TRUE, cores = 3)

  save(
    flash_model_ztnbinom,
    flash_model_ztSICHEL,
    file = "FlashAustriaModel.rda"
  )
}

## Compute quantile residuals.
resids <- c(
	"ztnbinom" = residuals(flash_model_ztnbinom, newdata = FlashAustriaTrain),
	"ztSICHEL" = residuals(flash_model_ztSICHEL, newdata = FlashAustriaTrain)
)

## Compare using worm plots.
plot(resids, which = "wp", main = "Worm plot",
  col = hcl.colors(2, "Dark 3"), cex = 0.5, ylim2 = c(-0.5, 0.5)) 

## Show coefficient paths of the boosting iterations of the best model.
pathplot(flash_model_ztSICHEL, which = "loglik.contrib", intercept = TRUE)

## Show some samples.
b2 <- flash_model_ztSICHEL
b2$samples[[1]] <- b2$samples[[1]][, c("mu.s.s(q_prof_PC1).b1", "mu.s.s(q_prof_PC1).b2")]
plot(b2, which = "samples", ask = FALSE)

## Plot selection of estimated effects.
plot(flash_model_ztSICHEL,
  term = c("s(q_prof_PC1)", "s(sqrt_cape)", "s(d2m)"),
  ask = FALSE, rug = TRUE, col.rug = "#39393919")

## Projection issues may happen, in this case try
sf::st_crs(FlashAustriaCase) <- 4326

## Compute estimated probabilities.
fit <- predict(flash_model_ztnbinom, newdata = FlashAustriaCase, type = "parameter")
fam <- family(flash_model_ztnbinom)
FlashAustriaCase$P10 <- 1 - fam$p(9, fit)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() + geom_sf(aes(fill = P10), data = FlashAustriaCase) +
  scale_fill_continuous_sequential("Oslo", rev = TRUE, p1 = .85) +
  geom_sf(data = world, col = "white", fill = NA) +
  coord_sf(xlim = c(7.95, 17), ylim = c(45.45, 50), expand = FALSE) +
  facet_wrap(~time, nrow = 2) + theme_minimal() +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))


# ---- Appendix A: Gaussian family object --------------------------------
Gauss_bamlss <- function(...) {
  f <- list(
    "family" = "mygauss",
    "names"  = c("mu", "sigma"),
    "links"  = c(mu = "identity", sigma = "log"),
    "d" = function(y, par, log = FALSE) {
      dnorm(y, mean = par$mu, sd = par$sigma, log = log)
    },
    "p" = function(y, par, ...) {
      pnorm(y, mean = par$mu, sd = par$sigma, ...)
    },
    "r" = function(n, par) {
      rnorm(n, mean = par$mu, sd = par$sigma)
    },
    "q" = function(p, par) {
      qnorm(p, mean = par$mu, sd = par$sigma)
    },
    "score" = list(
      mu = function(y, par, ...) {
        drop((y - par$mu) / (par$sigma^2))
      },
      sigma = function(y, par, ...) {
        drop(-1 + (y - par$mu)^2 / (par$sigma^2))
      }
    ),
    "hess" = list(
      mu = function(y, par, ...) {
        drop(1 / (par$sigma^2))
      },
      sigma = function(y, par, ...) { 
        rep(2, length(y))
      }
    )
  )
  class(f) <- "family.bamlss"
  return(f)
}



# ---- Appendix B: Special model terms -----------------------------------
## Smooth construct method.
smooth.construct.gc.smooth.spec <- function(object, data, knots) 
{
  object$X <- matrix(as.numeric(data[[object$term]]), ncol = 1)
  center <- if(!is.null(object$xt$center)) {
    object$xt$center
  } else TRUE
  object$by.done <- TRUE
  if(object$by != "NA")
    stop("by variables not supported!")

  ## Begin special elements to be used with bfit() and GMCMC().
  object$fit.fun <- function(X, b, ...) {
    f <- b[1] * exp(-b[2] * exp(-b[3] * drop(X)))
    if(center)
      f <- f - mean(f)
    f
  }
  object$update <- bfit_optim
  object$propose <- GMCMC_slice
  object$prior <- function(b) { sum(dnorm(b, sd = 1000, log = TRUE)) }
  object$fixed <- TRUE
  object$state$parameters <- c("b1" = 0, "b2" = 0.5, "b3" = 0.1)
  object$state$fitted.values <- rep(0, length(object$X))
  object$state$edf <- 3
  object$special.npar <- 3 ## Important!
  ## End special elements.

  ## Important, This is a special smooth constructor!
  class(object) <- c("gc.smooth", "no.mgcv", "special")

  object
}

## Predict matrix method.
Predict.matrix.gc.smooth <- function(object, data, knots) 
{
  X <- matrix(as.numeric(data[[object$term]]), ncol = 1)
  X
}

## Example using the growth curve smooth constructor
## using simulated data.
set.seed(111)

d <- data.frame("time" = 1:30)
d$y <- 2 + 1 / (1 + exp(0.5 * (15 - d$time))) +
  rnorm(30, sd = exp(-3 + 2 * cos(d$time/30 * 6 - 3)))

f <- list(
  y ~ s2(time, bs = "gc"),
  sigma ~ s(time)
)

b <- bamlss(f, data = d, optimizer = bfit, sampler = GMCMC)

## Plot estimated effects.
par(mfrow = c(1, 2))
p <- predict(b, model = "mu", FUN = c95)
plot(y ~ time, data = d, main = expression(mu))
plot2d(p ~ time, data = d, add = TRUE, axes = FALSE,
  fill.select = c(0, 1, 0, 1))
points(d$time, d$y)
plot(b, model = "sigma", spar = FALSE, main = expression(log(sigma)))


# ---- Appendix C: Model fitting engines for linear regression -----------
## Linear model family object.
lm_bamlss <- function(...) {
  f <- list(
    "family" = "LM",
    "names" = "mu",
    "links" = "identity",
    "d" = function(y, par, log = FALSE) {
      sigma <- sqrt(sum((y - par$mu)^2) / (length(y) - .lm_bamlss.p))
      dnorm(y, mean = par$mu, sd = sigma, log = log)
    },
    "p" = function(y, par, ...) {
      sigma <- sqrt(sum((y - par$mu)^2) / (length(y) - .lm_bamlss.p))
      pnorm(y, mean = par$mu, sd = sigma, ...)
    }
  )
  class(f) <- "family.bamlss"
  return(f)
}

## Simulate some data.
d <- GAMart()

## Test family by creating a bamlss frame.
bf <- bamlss.frame(num ~ x1 + x2, data = d, family = "lm")
print(bf)

## The linear model design matrix.
head(bf$x$mu$model.matrix)

## The response.
head(bf$y)

## Set up an optimizer function that can be used with bamlss().
lm.opt <- function(x, y, ...)
{
  ## Only univariate response.
  y <- y[[1L]]

  ## For illustration this is easier to read.
  X <- x$mu$model.matrix

  ## Estimate model parameters.
  par <- drop(chol2inv(chol(crossprod(X))) %*% crossprod(X, y))

  ## Set parameter names.
  names(par) <- paste0("mu.p.", colnames(X))

  ## Return estimated parameters and fitted values.
  rval <- list(
    "parameters" = par,
    "fitted.values" = drop(X %*% par),
    "edf" = length(par),
    "sigma" = drop(sqrt(crossprod(y - X %*% par) / (length(y) - ncol(X))))
  )

  ## Set edf within .GlobalEnv for the
  ## loglik() function in the lm_bamlss() family.
  .lm_bamlss.p <<- length(par)

  return(rval)
}

## Test the optimizer on the simulated data.
## Note, no MCMC sampling.
f <- num ~ x1 + poly(x2, 5) + poly(x3, 5)
b <- bamlss(f, data = d, family = "lm", optimizer = lm.opt, sampler = FALSE)
summary(b)
nd <- data.frame("x2" = seq(0, 1, length = 100))
nd$p <- predict(b, newdata = nd, term = "x2")

## Plot estimetd effect of x2.
plot2d(p ~ x2, data = nd)

## Set up a sampling engine for linear models.
lm.mcmc <- function(x, y, start = NULL,
  n.iter = 12000, burnin = 2000, thin = 10,
  m = 0, M = 1e+05,
  a = 1, b = 1e-05,
  verbose = TRUE, ...)
{
  ## How many samples are saved?
  itrthin <- seq.int(burnin, n.iter, by = thin)
  nsaves <- length(itrthin)

  ## Only univariate response.
  y <- y[[1L]]

  ## For illustration this is easier to read.
  X <- x$mu$model.matrix

  ## Again, set edf within .GlobalEnv for the
  ## loglik() function in the lm_bamlss() family.
  .lm_bamlss.p <<- ncol(X)

  ## Number of observations and parameters.
  n <- length(y)
  p <- ncol(X)

  ## Matrix saving the samples.
  samples <- matrix(0, nsaves, p + 1L)

  ## Stick to the naming convention.
  pn <- paste0("mu.p.", colnames(X))
  colnames(samples) <- c(
    pn,      ## Regression coefficients and
    "sigma"  ## variance samples.
  )

  ## Setup coefficient vector,
  ## again, use correct names.
  beta <- rep(0, p)
  names(beta) <- pn
  sigma <- sd(y)

  ## Check for starting values obtained,
  ## e.g., from lm.opt() from above.
  if(!is.null(start)) {
    sn <- names(start)
    for(j in names(beta)) {
      if(j %in% sn)
        beta[j] <- start[j]
    }
  }

  ## Process prior information.
  m <- rep(m, length.out = p)
  if(length(M) < 2)
    M <- rep(M, length.out = p)
  if(!is.matrix(M))
    M <- diag(M)
  Mi <- solve(M)

  ## Precompute cross products.
  XX <- crossprod(X)
  Xy <- crossprod(X, y)

  ## Inverse gamma parameter.
  a <- a + n / 2 + p / 2

  ## Start sampling.
  ii <- 1
  for(i in 1:n.iter) {
    ## Sampling sigma
    b2 <- b + 1 / 2 * t(y - X %*% beta) %*% (y - X %*% beta) +
      1 / 2 * t(beta - m) %*% Mi %*% (beta - m)
    sigma2 <- sqrt(1 / rgamma(1, a, b2))

    ## Sampling beta.
    sigma2i <- 1 / sigma2
    Sigma <- chol2inv(chol(sigma2i * XX + sigma2i * Mi))
    mu <- Sigma %*% (sigma2i * Xy + sigma2i * Mi %*% m)
    beta <- MASS::mvrnorm(1, mu, Sigma)
      
    if(i %in% itrthin) {
      samples[ii, pn] <- beta
      samples[ii, "sigma"] <- sqrt(sigma2)
      ii <- ii + 1
    }
    if(verbose) {
      if(i %% 1000 == 0)
        cat("iteration:", i, "\n")
    }
  }

  ## Convert to "mcmc" object.
  samples <- as.mcmc(samples)

  return(samples)
}

## Test both engines on the simulated data.
b <- bamlss(f, data = d, family = "lm", optimizer = lm.opt, sampler = lm.mcmc)

## Show the summary.
summary(b)

## Predict for all terms including 95% credible intervals
nd$x1 <- nd$x3 <- seq(0, 1, length = 100)
for(j in c("x1", "x2", "x3")) {
  nd[[paste0("p.", j)]] <- predict(b, newdata = nd, term = j,
    FUN = c95, intercept = FALSE)
}

## Plot estimated effects.
par(mfrow = c(1, 3))
plot2d(p.x1 ~ x1, data = nd, fill.select = c(0, 1, 0, 1), lty = c(2, 1, 2))
plot2d(p.x2 ~ x2, data = nd, fill.select = c(0, 1, 0, 1), lty = c(2, 1, 2))
plot2d(p.x3 ~ x3, data = nd, fill.select = c(0, 1, 0, 1), lty = c(2, 1, 2))

