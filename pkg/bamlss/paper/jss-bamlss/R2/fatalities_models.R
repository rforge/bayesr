## Required packages.
library("bamlss")
library("gamlss.dist")
library("parallel")

## Load the data.
data("fatalities", package = "bamlss")

## Full model formula.
f <- list(
  num   ~ s(week, bs = "cc"),
  sigma ~ s(week, bs = "cc"),
  nu    ~ s(week, bs = "cc"),
  tau   ~ s(week, bs = "cc")
)

## Family list.
families <- gamlss_distributions(type = "continuous")

## Setup function to run in parallel.
parallel_fun <- function(j) {
  cat("family", j, "\n")

  b1 <- try(bamlss(num ~ 1, data = fatalities,
    family = families[[j]](mu.link = "log"),
    optimizer = opt_boost, init = FALSE, plot = FALSE),
    silent = TRUE)

  b2 <- try(bamlss(f, data = fatalities,
    family = families[[j]](mu.link = "log"),
    optimizer = opt_boost, init = TRUE, plot = FALSE),
    silent = TRUE)

  rval <- list()
  rval$distribution <- j

  if(!inherits(b1, "try-error")) {
    par <- predict(b1, type = "parameter", drop = FALSE)
    dic <- DIC(b1)
    dnum <- family(b1)$d(fatalities$num, par)
    cat(".. .. b1: DIC =", dic$DIC, "pd =", dic$pd, "\n")
    rval$dic1 <- dic
    rval$dnum <- dnum
  } else {
    writeLines(b1)
  }

  if(!inherits(b2, "try-error")) {
    dic <- DIC(b2)
    cat(".. .. b2: DIC =", dic$DIC, "pd =", dic$pd, "\n")
    rval$dic2 <- dic
  } else {
    writeLines(b2)
  }

  return(rval)
}

## Estimate models.
res <- mclapply(names(families)[1:2], parallel_fun, mc.cores = 2)##length(families))
save(res, file = "fatalities_models.rda")

