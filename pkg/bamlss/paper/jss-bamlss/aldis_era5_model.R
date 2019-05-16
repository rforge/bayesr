## --- setup ---
set.seed(111)

## --- libraries ---
library("bamlss")

## --- data ---
load("ALDIS_ERA5_subset.rda")

## --- helper function to build formula ---
make_formula <- function(x, y, prefix = "s(", suffix = ", bs = 'ps')",
                         env = parent.frame()) {
  x <- x[x != y]
  x <- paste0(prefix, x, suffix)
  x <- paste0(x, collapse = " + ")
  f <- as.formula(paste("~", x))
  f <- list(update(f, substitute(a ~ ., list(a = as.name(y)))), f)
  environment(f) <- env
  f
}

## --- stability selection ---
f <- make_formula(names(d_train), "counts")
sel <- stabsel(f, data = d_train, family = "ztnbinom", q = 16, B = 100)

## --- fit selected model ---
newf <- formula(sel)
fam <- family(sel)

b <- bamlss(newf, data = d_train,             # standard interface
  family = fam, binning = TRUE,               # general arguments
  optimizer = boost, maxit = 1000,            # boosting arguments
  thin = 5, burnin = 1000, n.iter = 6000      # sampler arguments
)

## --- out-of-sample prediciton ---
fit <- predict(b, newdata = d_eval, type = "parameter")
d_eval <- cbind(d_eval, as.data.frame(fit))

## --- save results ---
save(b, d_eval, sel, file = "full_model.rda")



