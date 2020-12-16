## --- setup ---
set.seed(111)

## --- libraries ---
library("bamlss")
library("gamlss.dist")
library("gamlss.tr")

## --- data ---
load("ALDIS_ERA5_subset.rda")

## --- generate zero truncated Sichel distribution ---
SItr <- trun(0, family = "SI", local = FALSE)

## --- helper function to build formula ---
make_formula <- function(x, y, prefix = "s(", suffix = ", bs = 'ps')",
                         env = parent.frame(), k = 3) {
  x <- x[x != y]
  x <- paste0(prefix, x, suffix)
  x <- paste0(x, collapse = " + ")
  f <- as.formula(paste("~", x))
  f <- c(list(update(f, substitute(a ~ ., list(a = as.name(y))))),
    rep(list(f), k - 1))
  environment(f) <- env
  f
}

## --- stability selection ---
f <- make_formula(names(d_train), "counts")

if(!file.exists("sel_SItr.rda")) {
  sel_SItr <- stabsel(f, data = d_train, family = SItr, q = 16,
    B = 100L, thr = 0.9, maxit =  400, fraction = 0.25, cores = 50)
  save(sel_SItr, file = "sel_SItr.rda")
} else {
  load("sel_SItr.rda")
}

if(!file.exists("sel_ztnbinom.rda")) {
  sel_ztnbinom <- stabsel(f, data = d_train, family = "ztnbinom", q = 16,
    B = 100L, thr = 0.9, maxit =  400, fraction = 0.25, cores = 50)
  save(sel_ztnbinom, file = "sel_ztnbinom.rda")
} else {
  load("sel_ztnbinom.rda")
}

## --- fit selected models ---
newf1 <- formula(sel_SItr)
newf2 <- formula(sel_ztnbinom)

if(!file.exists("b1.rda")) {
  b1 <- bamlss(newf1, data = d_train,        # standard interface
    family = SItr, binning = TRUE,           # general arguments
    optimizer = opt_boost, maxit = 1000,         # boosting arguments
    thin = 5, burnin = 1000, n.iter = 6000,  # sampler arguments
    init = FALSE, light = TRUE
  )

  save(b1, file = "b1.rda")
} else {
  load("b1.rda")
}

if(!file.exists("b2.rda")) {
  b2 <- bamlss(newf2, data = d_train,
    family = "ztnbinom", binning = TRUE,
    optimizer = opt_boost, maxit = 1000,
    thin = 5, burnin = 1000, n.iter = 6000,
    init = FALSE, light = TRUE
  )

  save(b2, file = "b2.rda")
} else {
  load("b2.rda")
}

## --- predictions ---
d_eval <- na.omit(d_eval)
d_train <- na.omit(d_train)

fit1 <- predict(b1, newdata = d_eval, type = "parameter")
fit2 <- predict(b2, newdata = d_eval, type = "parameter")

family(b1)$loglik(d_eval$counts, fit1)
family(b2)$loglik(d_eval$counts, fit2)

e1t <- residuals(b1, newdata = d_train)
e2t <- residuals(b2, newdata = d_train)

e1e <- residuals(b1, newdata = d_eval)
e2e <- residuals(b2, newdata = d_eval)

par(mfrow = c(1, 2))
plot(c("SItr" = e1t, "ztnbinom" = e2t), which = "wp", main = "Training", pos = "top")
plot(c("SItr" = e1e, "ztnbinom" = e2e), which = "wp", main = "Testing", pos = "top")

## b1
logLik -36644.9 eps 0.0040 iteration 1000 qsel 26
 elapsed time: 116.92min
|********************| 100%  0.00sec 2880.81min

## b2
logLik -36511.1 eps 0.0006 iteration 1000 qsel 24
 elapsed time:  5.66min
|********************| 100%  0.00sec 83.89min



