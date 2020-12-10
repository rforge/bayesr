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
  sel_SItr <- stabsel(f, data = d_train, family = SItr, q = 16, B = 100, thr = .1)
  save(sel_SItr, file = "sel_SItr.rda")
} else {
  load("sel_SItr.rda")
}

if(!file.exists("sel_ztnbinom.rda")) {
  sel_ztnbinom <- stabsel(f, data = d_train, family = "ztnbinom", q = 16, B = 100, thr = .1)
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
    optimizer = boost, maxit = 1000,         # boosting arguments
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
    optimizer = boost, maxit = 1000,
    thin = 5, burnin = 1000, n.iter = 6000,
    init = FALSE, light = TRUE
  )

  save(b2, file = "b2.rda")
} else {
  load("b2.rda")
}

## --- out-of-sample prediciton ---
d_eval <- na.omit(d_eval)
fit1 <- predict(b1, newdata = d_eval, type = "parameter")
fit2 <- predict(b2, newdata = d_eval, type = "parameter")

d_eval <- cbind(d_eval,
  "SItr" = fit1,
  "ztnbinom" = fit2
)

## --- save results ---
save(b1, b2, d_eval, sel_SItr, sel_ztnbinom, file = "full_model.rda")

fit1 <- predict(b1, newdata = d_eval, type = "parameter")
fit2 <- predict(b2, newdata = d_eval, type = "parameter")

family(b1)$loglik(d_eval$counts, fit1)
family(b2)$loglik(d_eval$counts, fit2)

e1 <- residuals(b1, newdata = d_eval)
e2 <- residuals(b2, newdata = d_eval)

plot(c(e1, e2), which = "wp")

