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

## Setup function to run in parallel.
parallel_fun <- function(j) {
  cat("family", j, "\n")

  fam <- get(j)

  b1 <- try(bamlss(num ~ 1, data = fatalities,
    family = fam(mu.link = "log"), optimizer = opt_boost, init = FALSE,
    plot = FALSE, maxit = 1000, n.iter = 12000, burnin = 2000, thin = 10),
    silent = TRUE)

  b2 <- try(bamlss(f, data = fatalities,
    family = fam(mu.link = "log"), optimizer = opt_boost, init = TRUE,
    plot = FALSE, maxit = 1000, n.iter = 12000, burnin = 2000, thin = 10),
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
    rval$b1 <- b1
  } else {
    writeLines(b1)
  }

  if(!inherits(b2, "try-error")) {
    dic <- DIC(b2)
    cat(".. .. b2: DIC =", dic$DIC, "pd =", dic$pd, "\n")
    rval$dic2 <- dic
    rval$b2 <- b2
  } else {
    writeLines(b2)
  }

  return(rval)
}

## Families.
families <- c("NO", "GA", "BCT", "JSU", "BCPE", "BCCG")

## Estimate models.
res <- mclapply(families, parallel_fun, mc.cores = length(families))

## Save results.
save(res, file = "fatalities_models.rda")

