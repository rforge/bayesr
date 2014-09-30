library("BayesR")
library("lattice")

dgp2d <- function(n, sd = 0.2, type = c("simple", "complicated"))
{
  if(type == "simple") {
    f <- simfun(type = "2d")
  } else {
    f <- function(x1, x2) {
      xf <- scale2(simfun(type = "pick")(x1) * simfun("sin")(x2), -0.5, 0.5)
      xf <- xf - mean(xf)
     xf
    }
  }
  dat <- data.frame("x1" = sort(runif(n, 0, 1)), "x2" = runif(n, 0, 1))
  dat$f <- with(dat, f(x1, x2))
  dat$y <- with(dat, 1.2 + f + rnorm(n, sd = sd))
  dat
}

sim2d <- function(nrep = 100,
  n = c(100, 500, 1000),
  type = c("simple", "complicated"),
  sd = c(0.1, 0.2, 0.5),
  k = 10, sp = 0, fx = TRUE,
  plot = FALSE)
{
  sce <- expand.grid("rep" = 1:nrep, "n" = n, "sd" = sd,
    "type" = type, "k" = k, "sp" = sp, "fx" = fx)
  sce <- sce[order(sce$rep), ]

  res <- NULL; ii <- 1
  cat("** replication = 1  i = 1", rep(" ", 5))
  for(i in 1:nrow(sce)) {
    cat("\r")
    cat(paste("** replication =", sce$rep[i]), "i =", ii, rep(" ", 10))
    d <- dgp2d(n = sce$n[i], sd = sce$sd[i], type = sce$type[i])
    b1 <- gam(y ~ s(x1,x2,k=sce$k[i]*2,fx=sce$fx[i]), data = d)
    f <- as.formula(paste('y ~ rs(s(x1,k=', sce$k[i], ',sp=', 1/sp,
      ',fx=', sce$fx[i], '), s(x2,k=', sce$k[i], ',sp=', 1/sp, ',fx=', sce$fx[i],
      '), link="inverse")', sep = ''))
    b2 <- bayesr(f, data = d, method = "MP2", n.samples = 0, verbose = FALSE)

    f1 <- drop(predict(b1, type = "terms"))
    f1 <- f1 - mean(f1)
    f2 <- predict(b2, model = 1, term = 1, intercept = FALSE)
    f2 <- f2 - mean(f2)

    if(plot) {
      par(mfrow = c(1, 3))
      d$p1 <- f1; d$p2 <- f2
      plot3d(p1 ~ x1 + x2, data = d, main = paste("gam", sce$type[i]), type = "mba")
      plot3d(p2 ~ x1 + x2, data = d, main = paste("rs", sce$type[i]), type = "mba")
      plot3d(f ~ x1 + x2, data = d, main = paste("truth", sce$type[i]), type = "mba")
    }

    mse1 <- mean((d$f - f1)^2)
    mse2 <- mean((d$f - f2)^2)
    ll1 <- as.numeric(logLik(b1)[1])
    ll2 <- as.numeric(logLik(b2)[1])

    tres <- data.frame("rep" = rep(i, 2), "n" = rep(sce$n[i], 2), "sd" = rep(sce$sd[i], 2),
      "logLik" = c(ll1, ll2), "mse" = c(mse1, mse2), "model" = c("gam", "rs"), 
      "k" = rep(sce$k[i], 2), "fx" = rep(sce$fx[i], 2), "type" = rep(sce$type[i], 2))
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

res2d <- sim2d()

bwplot(log(sqrt(mse)) ~ model | sd + n, data = res2d)

