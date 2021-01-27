## Required packages.
library("bamlss")
library("gamlss.dist")
library("parallel")

## Load the data.
data("fatalities", package = "bamlss")

## Subset data 2000-2019.
d19 <- subset(fatalities, year < 2020)

## Full model formula.
f <- list(
  num   ~ s(week, bs = "cc", k = 20),
  sigma ~ s(week, bs = "cc", k = 20),
  nu    ~ s(week, bs = "cc", k = 20),
  tau   ~ s(week, bs = "cc", k = 20)
)

## Setup function to run in parallel.
parallel_fun <- function(j) {
  cat("family", j, "\n")

  fam <- get(j)

  set.seed(123)
  b1 <- bamlss(num ~ 1, data = d19, family = fam(mu.link = "log"),
    n.iter = 12000, burnin = 2000, thin = 10)

  set.seed(456)
  b2 <- bamlss(f, data = d19, family = fam(mu.link = "log"),
    n.iter = 12000, burnin = 2000, thin = 10)

  set.seed(789)
  k <- 10
  folds <- rep(1:k, length.out = nrow(d19))
  crps <- NULL
  for(i in 1:k) {
    df <- subset(d19, folds != i)
    de <- subset(d19, folds == i)
    b3 <- bamlss(f, data = df, family = fam(mu.link = "log"),
      n.iter = 12000, burnin = 2000, thin = 10)
    crps <- c(crps, CRPS(b3, newdata = de, FUN = identity))
  }

  rval <- list()
  rval$distribution <- j

  par <- predict(b1, type = "parameter", drop = FALSE)
  dic <- DIC(b1)
  waic <- WAIC(b1)
  dnum <- family(b1)$d(d19$num, par)
  cat(".. .. b1: DIC =", dic$DIC, "pd =", dic$pd, "WAIC =", waic$WAIC1, "\n")
  rval$dic1 <- dic
  rval$waic1 <- waic$WAIC1
  rval$dnum <- dnum
  rval$b1 <- b1

  dic <- DIC(b2)
  waic <- WAIC(b2)
  cat(".. .. b2: DIC =", dic$DIC, "pd =", dic$pd, "WAIC =", waic$WAIC1, "\n")
  rval$dic2 <- dic
  rval$waic2 <- waic$WAIC1
  rval$b2 <- b2
  rval$crps <- mean(crps)

  return(rval)
}

## Families.
families <- c("NO", "GA", "GG", "BCT", "BCPE", "BCCG")

## Estimate models.
res <- mclapply(families, parallel_fun, mc.cores = length(families))

## Save results.
save(res, file = "fatalities_models.rda")

