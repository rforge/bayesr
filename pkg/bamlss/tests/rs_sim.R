## Required packages.
library("bamlss")
library("lattice")

## Data generating process.
dgp <- function(n = 500, sd = 0.2, type = c("simple2d", "complicated2d"))
{
  type <- type[1]
  if(grepl("2d", type)) {
    if(type == "simple2d") {
      f <- simfun(type = "2d")
    } else {
      f <- function(x1, x2) {
        xf <- scale2(simfun(type = "pick")(x1) * simfun("sin")(x2), -0.5, 0.5)
        xf <- xf - mean(xf)
       xf
      }
    }
  } else {
    f <- function(x, ...) {
      simfun(type)(x)
    }
  }
  dat <- data.frame("x1" = sort(runif(n, 0, 1)), "x2" = runif(n, 0, 1))
  dat$ftrue <- with(dat, f(x1, x2))
  dat$y <- with(dat, 1.2 + ftrue + rnorm(n, sd = sd))
  dat
}

## Simulation and evaluation.
sim <- function(nrep = 100,
  n = c(100, 500, 1000),
  type = c("simple2d", "complicated2d", "sinus", "complicated"),
  sd = c(0.1, 0.2, 0.5),
  k = 10, sp = 0, fx = TRUE,
  plot = FALSE)
{
  sce <- expand.grid("rep" = 1:nrep, "n" = n, "sd" = sd,
    "type" = type, "k" = k, "sp" = sp, "fx" = fx, stringsAsFactors = FALSE)
  sce <- sce[order(sce$rep), ]

  res <- NULL; ii <- 1
  cat("** replication = 1  i = 1", "of", nrow(sce), rep(" ", 5))
  for(i in 1:nrow(sce)) {
    cat("\r")
    cat(paste("** replication =", sce$rep[i]), "i =", ii, "of", nrow(sce), rep(" ", 5))
    d <- dgp(n = sce$n[i], sd = sce$sd[i], type = sce$type[i])

    if(grepl("2d", sce$type[i])) {
      b1 <- gam(y ~ s(x1,x2,k=sce$k[i]*2,fx=sce$fx[i]), data = d)
      f <- as.formula(paste('y ~ rs(s(x1,k=', sce$k[i], ',sp=', 1/sp,
        ',fx=', sce$fx[i], '), s(x2,k=', sce$k[i], ',sp=', 1/sp, ',fx=', sce$fx[i],
        '), link="inverse")', sep = ''))
      b2 <- bamlss(f, data = d, method = "backfitting", n.samples = 0, verbose = FALSE)
    } else {
      b1 <- gam(y ~ s(x1,k=sce$k[i]*2,fx=sce$fx[i],bs="ps"), data = d)
      f <- as.formula(paste('y ~ rs(s(x1,k=', sce$k[i], ',sp=', 1/sp,
        ',fx=', sce$fx[i], ',bs="ps"), link="log")', sep = ''))
      b2 <- bamlss(f, data = d, method = "backfitting", n.samples = 0, verbose = FALSE)
    }

    f1 <- drop(predict(b1, type = "terms"))
    f1 <- f1 - mean(f1)
    f2 <- predict(b2, model = 1, term = 1, intercept = FALSE)
    f2 <- f2 - mean(f2)

    if(plot) {
      d$f1 <- f1; d$f2 <- f2
      if(grepl("2d", sce$type[i])) {
        par(mfrow = c(1, 3))
        plot3d(f1 ~ x1 + x2, data = d, main = paste("gam", sce$type[i]), type = "mba")
        plot3d(f2 ~ x1 + x2, data = d, main = paste("rs", sce$type[i]), type = "mba")
        plot3d(ftrue ~ x1 + x2, data = d, main = paste("truth", sce$type[i]), type = "mba")
      } else {
        par(mfrow = c(1, 1))
        plot2d(f1 ~ x1, data = d, main = paste("type", sce$type[i]), col.lines = 1, rug = FALSE)
        plot2d(f2 ~ x1, data = d, add = TRUE, col.lines = 3, rug = FALSE)
        plot2d(ftrue ~ x1, data = d, add = TRUE, col.lines = 2, rug = FALSE)
        legend("topleft", c("truth", "gam", "rs"), col = c(2, 1, 3), lwd = 2)
      }
    }

    rmse1 <- sqrt(mean(residuals(b1)^2))
    rmse2 <- sqrt(mean(residuals(b2, type = "ordinary")^2))
    bias1 <- mean((d$ftrue - f1)^2)
    bias2 <- mean((d$ftrue - f2)^2)
    ll1 <- as.numeric(logLik(b1)[1])
    ll2 <- as.numeric(logLik(b2)[1])

    tres <- data.frame("rep" = rep(i, 2), "n" = rep(sce$n[i], 2), "sd" = rep(sce$sd[i], 2),
      "logLik" = c(ll1, ll2), "bias" = c(bias1, bias2), "rmse" = c(rmse1, rmse2),
      "model" = c("gam", "rs"), "k" = rep(sce$k[i], 2), "fx" = rep(sce$fx[i], 2),
      "type" = rep(sce$type[i], 2))
    res <- rbind(res, tres)

    ii <- ii + 1
  }

  cat("\n")

  res$sd <- as.factor(res$sd)
  res$n <- as.factor(res$n)
  res$fx <- as.factor(res$fx)
  res$k <- as.factor(res$k)
  res$type <- as.factor(res$type)

  return(res)
}


if(FALSE) {
  rsim <- sim(nrep = 100)
  save(rsim, file = "rsim.rda")
  bwplot(log(sqrt(bias)) ~ model | sd + n, data = rsim)
}

