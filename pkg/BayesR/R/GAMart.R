########################################
## (1) Functions used for simulation. ##
########################################
simfun <- function(type = "sinus")
{
  ## Known function types.
  known_types <- c("linear", "quadratic", "unimodal", "double", "sinus",
    "cosinus", "pick", "complicated", "const", "spatial")
  if(is.character(type))
    type <- match.arg(type, known_types)

  f <- switch(type,
    "linear" = function(x, min = 0, max = 1) {
      x <- (x - min) / (max - min)
      -4.5 * x
    },
    "quadratic" = function(x, min = 0, max = 1) {
      x <- (x - min) / (max - min)
      3.5 * (x - 0.5)^2
    },
    "unimodal" = function(x, min = 0, max = 1) {
      x <- (x - min) / (max - min)
      120 * x * exp(-10 * x)
    },
    "double" = function(x, min = 0, max = 1) {
      x <- (x - min) / (max - min)
      1.3 * (120 * x * exp(-10 * x) + 2.75 * x^2)
    },
    "sinus" = function(x, min = 0, max = 1) {
      x <- (x - min) / (max - min)
      sin(2 * pi * x)
    },
    "cosinus" = function(x, min = 0, max = 1) {
      x <- (x - min) / (max - min)
      cos(2 * pi * x)
    },
    "pick" = function(x, min = 0, max = 1) { 
      x <- (x - min) / (max - min)
      sin(2 * (4 * x - 2)) + 2 * exp(-16^2 * (x - 0.5)^2)
    },
    "complicated" = function(x, min = 0, max = 1) {
      x <- (x - min) / (max - min)
      exp(-400 * (x - 0.6)^2) + 5 * exp(-500 * (x - 0.75)^2) / 3 + 2 * exp(-500 * (x - 0.9)^2)
    },
    "const" = function(..., const = 1.2) {
      const
    },
    "spatial" = function(id, f1 = sin, f2 = function(x) sin(x * pi)) {
      n <- ceiling(sqrt(length(id)))
      co <- expand.grid("long" = seq(0, 1, length = n), "lat" = seq(0, 1, length = n))
      f <- f1(co[, 1]) * f2(co[, 2])
      f <- data.frame("long" = co[, 1], "lat" = co[, 2], "f" = f)
      f[seq_along(id), ]
    }
  )
  if(!is.character(type))
    type <- known_types[type]
  attr(f, "type") <- type

  f
}

## Function for scaling.
scale2 <- function(x, lower = -1.5, upper = 1.5)
{
  x <- if(length(unique(x)) > 1) {
    (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE)) * (upper - lower) + lower
  } else x
  x
}


#####################################################
## (2) Hierarchical predictor generating process.  ##
#####################################################
dgp_eta <- function(nobs = c(500, 100, 10), const = 1.2,
  type = list(c("unimodal", "linear"), c("sinus", "spatial"), "const"),
  scale = list(c(1, 1), c(1, 1)), sd = rep(0.3, length(nobs) - 1),
  nx = 100, min = 0, max = 1, seed0 = 1)
{
  ## Start seed, only needed for random effects.
  if(is.null(seed0))
    seed0 <- round(runif(1L) * .Machine$integer.max)

  ## Other seeds
  seed1 <- round(runif(length(nobs)) * .Machine$integer.max)

  ## Number of hierarchical levels.
  nlevels <- length(nobs)

  ## Create covariates and predictor for each level.
  d <- vector(mode = "list", length = nlevels)
  maps <- list()
  for(j in nlevels:1) {
    fun <- xd <- list()
    for(i in seq_along(type[[j]])) {
      if(any(duplicated(type[[j]]))) stop("duplicated function specification!")
      FUN <- simfun(type[[j]][[i]])
      mt <- attr(FUN, "type")
      if(i < 2 & j > 1) {
        xd[[id <- paste("id", j, sep = "")]] <- factor(1:nobs[j])
        set.seed(seed1[j])
        fun[[paste("re", j, sep = "")]] <- rnorm(nobs[j], sd = sd[j - 1])
      }
      if(mt == "spatial") {
        e <- FUN(1:nobs[j])
        xd[[paste("long", j, sep = "")]] <- e$long
        xd[[paste("lat", j, sep = "")]] <- e$lat
        fun[[paste("sp", j, sep = "")]] <- scale2(e$f) * if(is.na(scale[[j]][i])) 1 else scale[[j]][i]
        m <- pixelmap(lat ~ long, data = e, at.xy = TRUE, size = 0.5, yscale = FALSE)
        names(m$map) <- as.character(xd[[id]])
        m$data$id <- as.integer(as.character(xd[[id]]))
        maps[[paste("map", j, sep = "")]] <- m
      } else {
        xn <- paste("x", j, i, sep = "")
        fn <- if(mt != "const") paste("f", j, i, sep = "") else "const"
        set.seed(seed0 + j + i)
        x <- sample(rep(seq(min, max, length = nx), length.out = nobs[j]))
        fun[[fn]] <- if(mt != "const") {
          scale2(FUN(x)) * if(is.na(scale[[j]][i])) 1 else scale[[j]][i]
        } else rep(FUN(const = const), length.out = nobs[j])
        if(mt != "const") xd[[xn]] <- x
      }
    }
    eta <- 0
    fn <- names(fun)
    for(i in seq_along(fun)) {
      if(fn[i] != "const")
        fun[[i]] <- fun[[i]] - mean(fun[[i]])
      eta <- eta + fun[[i]]
    }
    fun[[paste("eta", j, sep = "")]] <- eta
    d[[j]] <- as.data.frame(c(xd, fun))
  }

  ## Wrap up to a final data.frame.
  d <- as.data.frame(d)

  ## Create final predictor.
  d$eta0 <- rowSums(d[, grepl("eta", names(d)), drop = FALSE])

  ## Assign maps.
  attr(d, "maps") <- maps

  return(d)
}


####################################
## (3) Data generating functions. ##
####################################
## Gaussian.
dgp_gaussian <- function(n = 500, mu = NULL, sigma = NULL, range.sigma = c(0.3, 1.2), ...)
{
  if(is.null(mu)) {
    mu <- list(nobs = n, const = 0.5,
      type = list(c("unimodal", "quadratic", "const")))
  }
  if(is.null(sigma)) {
    sigma <- list(nobs = n, const = 0.01,
      type = list(c("pick", "const")))
  }

  mu <- do.call("dgp_eta", mu)
  sigma <- do.call("dgp_eta", sigma)

  d <- data.frame("mu" = mu, "sigma" = sigma)
  sd <- (scale2(sigma$eta0, range.sigma[1], range.sigma[2]))
  d$y <- rnorm(nrow(mu), mean = mu$eta0, sd = sd)

  d
}

if(FALSE) {
  d <- dgp_gaussian()
  b <- bayesr(y ~ sx(mu.x11) + sx(mu.x12),y ~ sx(sigma.x11), data = d, engine = "BayesX", verbose = TRUE)
  d$p <- predict(b, model = "mu", term = c("x11", "x12"))
}

## Gaussian2.
dgp_gaussian2 <- function(n = 500, mu = NULL, sigma2 = NULL, range.sigma2 = c(0.3, 1.2), ...)
{
  if(is.null(mu)) {
    mu <- list(nobs = n, const = 0.5,
      type = list(c("unimodal", "quadratic", "const")))
  }
  if(is.null(sigma2)) {
    sigma2 <- list(nobs = n, const = 0.01,
      type = list(c("pick", "const")))
  }

  mu <- do.call("dgp_eta", mu)
  sigma2 <- do.call("dgp_eta", sigma2)

  d <- data.frame("mu" = mu, "sigma2" = sigma2)
  sd <- sqrt(scale2(sigma2$eta0, range.sigma2[1], range.sigma2[2]))
  d$y <- rnorm(nrow(mu), mean = mu$eta0, sd = sd)

  d
}

if(FALSE) {
  d <- dgp_gaussian()
  b <- bayesr(y ~ sx(mu.x11) + sx(mu.x12), y ~ sx(sigma2.x11), data = d, engine = "BayesX", verbose = TRUE)
  d$p <- predict(b, model = "mu", term = c("x11", "x12"))
}

## truncated gaussian2.
dgp_truncgaussian <- function(n = 500, mu = NULL, sigma = NULL, ...)
{
  if(is.null(mu)) {
    mu <- list(nobs = n, const = 0.5,
      type = list(c("unimodal", "quadratic", "const")))
  }
  if(is.null(sigma)) {
    sigma <- list(nobs = n, const = 0.01,
      type = list(c("const")))
  }

  mu <- do.call("dgp_eta", mu)
  sigma <- do.call("dgp_eta", sigma)

  d <- data.frame("mu" = mu, "sigma" = sigma)
  sd <- exp(sigma$eta0)
  u <- runif(n)
  d$y <- qnorm(0.5 + 0.5 * u) * sd + mu$eta0

  d
}

if(FALSE) {
  d <- dgp_truncgaussian()
  b <- bayesr(y ~ sx(mu.x11) + sx(mu.x12), y ~ 1, data = d, engine = "BayesX", verbose = TRUE, family = truncgaussian)
  d$p <- predict(b, model = "mu", term = c("x11", "x12"))
  plot(b, which = 3:6)
}

## truncated gaussian2.
dgp_truncgaussian2 <- function(n = 500, mu = NULL, sigma2 = NULL, ...)
{
  if(is.null(mu)) {
    mu <- list(nobs = n, const = 0.5,
      type = list(c("unimodal", "quadratic", "const")))
  }
  if(is.null(sigma2)) {
    sigma2 <- list(nobs = n, const = 0.01,
      type = list(c("const")))
  }

  mu <- do.call("dgp_eta", mu)
  sigma2 <- do.call("dgp_eta", sigma2)

  d <- data.frame("mu" = mu, "sigma2" = sigma2)
  sd <- sqrt(exp(sigma2$eta0))
  u <- runif(n)
  d$y <- qnorm(0.5 + 0.5 * u) * sd + mu$eta0

  d
}

if(FALSE) {
  d <- dgp_truncgaussian2()
  b <- bayesr(y ~ sx(mu.x11) + sx(mu.x12), y ~ 1, data = d, engine = "BayesX", verbose = TRUE, family = truncgaussian)
  d$p <- predict(b, model = "mu", term = c("x11", "x12"))
  plot(b, which = 3:6)
}


## Gamma.
dgp_gamma <- function(n = 500, mu = NULL, sigma = NULL, ...)
{
  if(is.null(mu)) {
    mu <- list(nobs = n, const = 0,
      type = list(c("unimodal", "quadratic", "const")))
  }
  if(is.null(sigma)) {
    sigma <- list(nobs = n, const = 0,
      type = list(c("linear", "sinus", "const")))
  }


  mu <- do.call("dgp_eta", mu)
  sigma <- do.call("dgp_eta", sigma) 
  m <- exp(mu$eta0)
  s <- exp(sigma$eta0)
  y <- rgamma(n, shape = s, scale = m/s)

  
  d <- cbind(y, "mu" = mu, "sigma" = sigma)
  d
}

if(FALSE) {
  d <- dgp_gamma()
  b <- bayesr(y ~ sx(mu.x11) + sx(mu.x12), y ~ sx(sigma.x11)+sx(sigma.x12), 
		data = d, family = gamma, engine = "BayesX", verbose = TRUE)
  d$pred_mu <- predict(b, model = "mu", term = c("x11", "x12"))
}

## Inverse Gaussian.
dgp_invgaussian <- function(n = 500, mu = NULL, sigma = NULL, ...)
{
  require(gamlss)
  if(is.null(mu)) {
    mu <- list(nobs = n, const = 0,
      type = list(c("unimodal", "quadratic", "const")))
  }
  if(is.null(sigma)) {
    sigma <- list(nobs = n, const = 0,
      type = list(c("linear", "sinus", "const")))
  }


  mu <- do.call("dgp_eta", mu)
  sigma <- do.call("dgp_eta", sigma) 
  m <- exp(mu$eta0)
  s <- exp(sigma$eta0)
  y <- rIG(n, mu = m, sigma = sqrt(s))

  
  d <- cbind(y, "mu" = mu, "sigma2" = sigma)
  d
}

if(FALSE) {
  d <- dgp_invgaussian()
  b <- bayesr(y ~ sx(mu.x11) + sx(mu.x12), y ~ sx(sigma2.x11)+sx(sigma2.x12), 
		data = d, family = invgaussian, engine = "BayesX", verbose = TRUE)
  d$pred_mu <- predict(b, model = "mu", term = c("x11", "x12"))
}

## Lognormal.
dgp_lognormal <- function(n = 500, mu = NULL, sigma = NULL, ...)
{
  if(is.null(mu)) {
    mu <- list(nobs = n, const = 0,
      type = list(c("complicated", "quadratic", "const")))
  }
  if(is.null(sigma)) {
    sigma <- list(nobs = n, const = 0,
      type = list(c("quadratic", "linear", "const")))
  }


  mu <- do.call("dgp_eta", mu)
  sigma <- do.call("dgp_eta", sigma) 
  m <- (mu$eta0)
  s <- exp(sigma$eta0)
  y <- rlnorm(n, meanlog = m, sdlog = (s))

  
  d <- cbind(y, "mu" = mu, "sigma" = sigma)
  d
}

if(FALSE) {
  d <- dgp_lognormal()
  b <- bayesr(y ~ sx(mu.x11) + sx(mu.x12), y ~ sx(sigma.x11)+sx(sigma.x12), 
		data = d, family = lognormal, engine = "BayesX")
  plot(b)
  summary(b)
  mu.f11.est <- predict(b, model = "mu", term = 1, what = "terms")
  score(b)
}


## Lognormal2.
dgp_lognormal2 <- function(n = 500, mu = NULL, sigma2 = NULL, ...)
{
  if(is.null(mu)) {
    mu <- list(nobs = n, const = 0,
      type = list(c("complicated", "quadratic", "const")))
  }
  if(is.null(sigma2)) {
    sigma2 <- list(nobs = n, const = 0,
      type = list(c("quadratic", "linear", "const")))
  }


  mu <- do.call("dgp_eta", mu)
  sigma2 <- do.call("dgp_eta", sigma2) 
  m <- (mu$eta0)
  s <- exp(sigma2$eta0)
  y <- rlnorm(n, meanlog = m, sdlog = sqrt(s))

  
  d <- cbind(y, "mu" = mu, "sigma2" = sigma2)
  d
}

if(FALSE) {
  d <- dgp_lognormal2()
  b <- bayesr(y ~ sx(mu.x11) + sx(mu.x12), y ~ sx(sigma2.x11)+sx(sigma2.x12), 
		data = d, family = lognormal2, engine = "BayesX")
  plot(b, which = 3:6)
  summary(b)
  score(b)
}


## weibull.
dgp_weibull <- function(n = 500, lambda = NULL, alpha = NULL, ...)
{
  if(is.null(lambda)) {
    lambda <- list(nobs = n, const = 0,
      type = list(c("unimodal", "quadratic", "const")))
  }
  if(is.null(alpha)) {
    alpha <- list(nobs = n, const = 0,
      type = list(c("const")))
  }


  lambda <- do.call("dgp_eta", lambda)
  alpha <- do.call("dgp_eta", alpha) 
  m <- exp(lambda$eta0)
  s <- exp(alpha$eta0)
  y <- rweibull(n, shape = s, scale = m)

  
  d <- cbind(y, "lambda" = lambda, "alpha" = alpha)
  d
}

if(FALSE) {
  d <- dgp_weibull()
  b <- bayesr(y ~ sx(lambda.x11) + sx(lambda.x12), y ~ 1, 
		data = d, family = weibull, engine = "BayesX", verbose = TRUE)
  plot(b)
}

## t.
dgp_t <- function(n = 1000, mu = NULL, sigma2 = NULL, df = NULL, ...)
{
  if(is.null(mu)) {
    mu <- list(nobs = n, const = 0.5,
      type = list(c("unimodal", "quadratic", "const")))
  }
  if(is.null(sigma2)) {
    sigma2 <- list(nobs = n, const = 0.01,
      type = list(c("pick", "const")))
  }
  if(is.null(df)) {
    df <- list(nobs = n, const = 0.01,
      type = list(c("linear", "const")))
  }

  mu <- do.call("dgp_eta", mu)
  sigma2 <- do.call("dgp_eta", sigma2)
  df <- do.call("dgp_eta", df)

  d <- data.frame("mu" = mu, "sigma2" = sigma2, "df" = df)
  sd <- sqrt(exp(sigma2$eta0))
  degf <- (exp(df$eta0))
  ytemp <- rt(n, df = degf)
  d$y <- mu$eta0 + sd*ytemp

  d
}

if(FALSE) {
  d <- dgp_t()
  b <- bayesr(y ~ sx(mu.x11) + sx(mu.x12), y ~ sx(sigma2.x11), y ~ sx(df.x11), data = d, family = t, engine = "BayesX", verbose = TRUE)
  plot(b)
}



## Dagum.
dgp_dagum <- function(n = 500, a = NULL, b = NULL, p = NULL, ...)
{
  require(VGAM)
  if(is.null(a)) {
    a <- list(nobs = n, const = 0,
      type = list(c("unimodal", "quadratic", "const")))
  }
  if(is.null(b)) {
    b <- list(nobs = n, const = 0,
      type = list(c("linear", "sinus", "const")))
  }
  if(is.null(p)) {
    p <- list(nobs = n, const = 0,
      type = list(c("const")))
  }

  a <- do.call("dgp_eta", a)
  b <- do.call("dgp_eta", b) 
  p <- do.call("dgp_eta", p)
  shape1.a <- exp(a$eta0)
  scale <- exp(b$eta0)
  shape2.p <- exp(p$eta0)
  y <- rdagum(n, shape1.a = shape1.a, scale = scale, shape2.p = shape2.p)

  
  d <- cbind(y, "a" = a, "b" = b, "p" = p)
  d
}

if(FALSE) {
  d <- dgp_dagum()
  b <- bayesr(list(y ~ sx(a.x11) + sx(a.x12), y~ sx(b.x11)+sx(b.x12), y~ 1),
		data = d, family = dagum, engine = "BayesX", verbose = TRUE)
  summary(b)
  plot(b)
}



## Beta.
dgp_beta <- function(n = 500, mu = NULL, sigma2 = NULL, ...)
{
	if(is.null(mu)) {
    mu <- list(nobs = n, const = -0.5,
      type = list(c("const")))
  }
  if(is.null(sigma2)) {
    sigma2 <- list(nobs = n, const = 0.01,
      type = list(c("const")))
  }
  
  mu <- do.call("dgp_eta", mu)
  sigma2 <- do.call("dgp_eta", sigma2)

  d <- data.frame("mu" = mu, "sigma2" = sigma2)
  b <- exp(sigma2$eta0)
  b <- b / (1 + b)
  a <- exp(mu$eta0)
  a <- a / (1 + a)
  shape1 <- a * (1 - b) / (b)
  shape2 <- (1 - a) * (1 - b) / (b)
  y <- rbeta(n, shape1 = shape1, shape2 = shape2)
  d <- cbind(y, "mu" = mu, "sigma2" = sigma2)
  d
}

if(FALSE) {
  d <- dgp_beta()
  b <- bayesr(y ~ 1, data = d, family = beta, engine = "BayesX", verbose = TRUE)
  summary(b)
  plot(b, which = 3:6)
}


## BCCG.
dgp_BCCG <- function(n = 1000, mu = NULL, sigma = NULL, nu = NULL, ...)
{
  require(gamlss)
  if(is.null(mu)) {
    mu <- list(nobs = n, const = 0.5,
      type = list(c("unimodal", "quadratic", "const")))
  }
  if(is.null(sigma)) {
    sigma <- list(nobs = n, const = 0.01,
      type = list(c("pick", "const")))
  }
  if(is.null(nu)) {
    nu <- list(nobs = n, const = -0.1,
      type = list(c("const")))
  }

  mu <- do.call("dgp_eta", mu)
  sigma <- do.call("dgp_eta", sigma)
  nu <- do.call("dgp_eta", nu)

  d <- data.frame("mu" = mu, "sigma" = sigma, "nu" = nu)
  mup <- (exp(mu$eta0))
  sd <- (exp(sigma$eta0))
  nup <- (nu$eta0)
  d$y <- rBCCG(n, mu = mup, sigma = sd, nu = nup)

  d
}

if(FALSE) {
  d <- dgp_BCCG()
  f <- list(
    y ~ sx(mu.x11) + sx(mu.x12),
    y ~ sx(sigma.x11),
	y ~ 1
  )
  b <- bayesr(f, data = d, family = BCCG, engine = "BayesX", verbose = TRUE, n.iter = 6000, burnin = 1000, thin = 5)
  plot(b)
}


## Multivariate normal.
dgp_mvn <- function(n = 1000, mu1 = NULL, mu2 = NULL, sigma1 = NULL, sigma2 = NULL, rho = NULL, ...)
{
  require("mvtnorm")

  if(is.null(mu1)) {
    mu1 <- list(nobs = n, const = 0,
      type = list(c("unimodal", "quadratic", "const")))
  }
  if(is.null(mu2)) {
    mu2 <- list(nobs = n, const = 0,
      type = list(c("linear", "sinus", "const")))
  }
  if(is.null(sigma1)) {
    sigma1 <- list(nobs = n, const = 0.01,
      type = list(c("pick", "const")))
  }
  if(is.null(sigma2)) {
    sigma2 <- list(nobs = n, const = 0.01,
      type = list(c("cosinus", "const")))
  }
  if(is.null(rho)) {
    rho <- list(nobs = n, const = 0.01,
      type = list(c("linear", "const")))
  }

  mu1 <- do.call("dgp_eta", mu1)
  mu2 <- do.call("dgp_eta", mu2)
  sigma1 <- do.call("dgp_eta", sigma1)
  sigma2 <- do.call("dgp_eta", sigma2)
  rho <- do.call("dgp_eta", rho)

  y <- NULL
  for(i in 1:nrow(mu1)) {
    m <- c(mu1$eta0[i], mu2$eta0[i])
    s1 <- exp(sigma1$eta0[i])
    s2 <- exp(sigma2$eta0[i])
    p <- rho$eta0[i] / sqrt(1 + rho$eta0[i]^2)
    sm <- matrix(c(s1^2, p * s1 * s2, p * s1 * s2, s2^2), 2, 2)
    y <- rbind(y, rmvnorm(1,mean = m, sigma = sm))
  }
  colnames(y) <- c("y1", "y2")
  y <- as.data.frame(y)
  
  d <- cbind(y, "mu1" = mu1, "mu2" = mu2, "sigma1" = sigma1, "sigma2" = sigma2, "rho" = rho)
  d
}

if(FALSE) {
  d <- dgp_mvn()

  f <- list(
    y1 ~ sx(mu1.x11) + sx(mu1.x12),
    y2 ~ sx(mu2.x11) + sx(mu2.x12),
    y1 ~ sx(sigma1.x11),
    y2 ~ sx(sigma2.x11),
    y1 ~ sx(rho.x11)
  )

  b <- bayesr(f, family = mvn, data = d, engine = "BayesX", verbose = TRUE)
  
  plot(b)
}


## zip.
dgp_zip <- function(n = 500, lambda = NULL, p = NULL, 
			range.lambda = c(-1,1), range.p = c(-1,1), ...)
{
  require("VGAM")

  if(is.null(lambda)) {
    lambda <- list(nobs = n, const = -0.5,
      type = list(c("unimodal", "quadratic", "const")))
  }
  if(is.null(p)) {
    p <- list(nobs = n, const = 0,
      type = list(c("linear", "sinus", "const")))
  }


  lambda <- do.call("dgp_eta", lambda)
  p <- do.call("dgp_eta", p) 
  ld <- exp(lambda$eta0)
  pi <- exp(p$eta0)
  pi <- pi/(1+pi)
  y <- rzipois(n, lambda = ld, pstr0 = pi)

  
  d <- cbind(y, "lambda" = lambda, "pi" = p)
  d
}

if(FALSE) {
  d <- dgp_zip()

  f <- list(
    y ~ sx(lambda.x11) + sx(lambda.x12),
    y ~ sx(pi.x11) + sx(pi.x12)
  )

  b <- bayesr(f, family = zip, data = d, engine = "BayesX", verbose = TRUE)
  lambda.f11.est <- predict(b, model = "lambda", term = 1, what = "terms")
  lambda.f11.est.2p5 <- predict(b, model = "lambda", term = 1, what = "terms", FUN = quantile, 0.025)
  lambda.f11.est.97p5 <- predict(b, model = "lambda", term = 1, what = "terms", FUN = quantile, 0.975)
  
  plot(d$lambda.x11[order(d$lambda.x11)], lambda.f11.est[order(d$lambda.x11)], type = "l", lty = 2)
  lines(d$lambda.x11[order(d$lambda.x11)], d$lambda.f11[order(d$lambda.x11)]-mean(d$lambda.f11[order(d$lambda.x11)]))
  lines(d$lambda.x11[order(d$lambda.x11)], lambda.f11.est.2p5[order(d$lambda.x11)], lty = 2)
  lines(d$lambda.x11[order(d$lambda.x11)], lambda.f11.est.97p5[order(d$lambda.x11)], lty = 2)
  
  
  plot(b)
  plot(b, which = 3:6)
  lambda <- predict(b, type = "parameter", model = "lambda")
  pi <- predict(b, type = "parameter", model = "pi")
  res <- qnorm(runif(length(d$y), min = pzipois(d$y - 1, lambda = lambda, pstr0 = pi), max = pzipois(d$y, lambda = lambda, pstr0 = pi)))
  qqnorm(res)
  
  
}


## poisson.
dgp_poisson <- function(n = 500, lambda = NULL, 
			range.lambda = c(-1,1), range.p = c(-1,1), ...)
{

  if(is.null(lambda)) {
    lambda <- list(nobs = n, const = -0.5,
      type = list(c("unimodal", "quadratic", "const")))
  }

  lambda <- do.call("dgp_eta", lambda)
  ld <- exp(lambda$eta0)
  y <- rpois(n, lambda = ld)

  
  d <- cbind(y, "lambda" = lambda)
  d
}

if(FALSE) {
  d <- dgp_poisson()

  f <- list(
    y ~ sx(lambda.x11) + sx(lambda.x12)
  )

  b <- bayesr(f, family = poisson, data = d, engine = "BayesX", verbose = TRUE)
  lambda.f11.est <- predict(b, model = "lambda", term = 1, what = "terms")
  lambda.f11.est.2p5 <- predict(b, model = "lambda", term = 1, what = "terms", FUN = quantile, 0.025)
  lambda.f11.est.97p5 <- predict(b, model = "lambda", term = 1, what = "terms", FUN = quantile, 0.975)
  
  plot(d$lambda.x11[order(d$lambda.x11)], lambda.f11.est[order(d$lambda.x11)], type = "l", lty = 2)
  lines(d$lambda.x11[order(d$lambda.x11)], d$lambda.f11[order(d$lambda.x11)]-mean(d$lambda.f11[order(d$lambda.x11)]))
  lines(d$lambda.x11[order(d$lambda.x11)], lambda.f11.est.2p5[order(d$lambda.x11)], lty = 2)
  lines(d$lambda.x11[order(d$lambda.x11)], lambda.f11.est.97p5[order(d$lambda.x11)], lty = 2)
  
  
  plot(b)
  
  
}


## negbin.
dgp_negbin <- function(n = 500, mu = NULL, delta = NULL, 
			range.mu = c(-1,1), range.delta = c(-1,1), ...)
{

  if(is.null(mu)) {
    mu <- list(nobs = n, const = -0.5,
      type = list(c("unimodal", "quadratic", "const")))
  }
  if(is.null(delta)) {
    delta <- list(nobs = n, const = 0,
      type = list(c("linear", "const")))
  }


  mu <- do.call("dgp_eta", mu)
  delta <- do.call("dgp_eta", delta) 
  ld <- exp(mu$eta0)
  d <- exp(delta$eta0)
  y <- rnbinom(n, mu = ld, size = d)

  
  d <- cbind(y, "mu" = mu, "delta" = delta)
  d
}

if(FALSE) {
  d <- dgp_negbin()

  f <- list(
    y ~ sx(mu.x11) + sx(mu.x12),
    y ~ (delta.x11) 
  )

  b <- bayesr(f, family = negbin, data = d, engine = "BayesX", verbose = TRUE)
  
  plot(b, which = 3:6)
  
  
}

## zinb.
dgp_zinb <- function(n = 500, mu = NULL, pi = NULL, delta = NULL, ...)
{
  require(VGAM)
  if(is.null(mu)) {
    mu <- list(nobs = n, const = -0.5,
      type = list(c("unimodal", "quadratic", "const")))
  }
  if(is.null(pi)) {
    pi <- list(nobs = n, const = 0.1,
      type = list(c("const")))
  }
  if(is.null(delta)) {
    delta <- list(nobs = n, const = 1,
      type = list(c("const")))
  }


  mu <- do.call("dgp_eta", mu)
  pi <- do.call("dgp_eta", pi) 
  delta <- do.call("dgp_eta", delta) 
  ld <- exp(mu$eta0)
  p <- exp(pi$eta0)
  p <- p / (1 + p)
  dl <- exp(delta$eta0)
  y <- rzinegbin(n, munb = ld, size = dl, pstr0 = p)

  
  d <- cbind(y, "mu" = mu, "pi" = pi, "delta" = delta)
  d
}

if(FALSE) {
  d <- dgp_zinb()

  f <- list(
    y ~ sx(mu.x11) + sx(mu.x12),
	y ~ 1,
    y ~  1
  )

  b <- bayesr(f, family = zinb, data = d, engine = "BayesX", verbose = TRUE)
  
  plot(b, which = 3:6)
  
  
}


## Multinomial.
dgp_multinomial <- function(nlevels = 4,
  nobs = c(1000, 100, 10),
  specs = list(
    type = list(c("unimodal", "linear"), c("sinus", "spatial"), "const"),
    scale = list(c(0.5, 1), c(0.7, 1)),
    sd = rep(0.3, length(nobs) - 1),
    const = 0
  ))
{
  ident <- FALSE
  if(!is.null(names(specs))) {
    specs <- list(specs)
    ident <- TRUE
  }
  specs <- rep(specs, length.out = nlevels - 1)

  d <- maps <- list()
  for(j in 2:nlevels) {
    specs[[j - 1]]$nobs <- nobs
    eta <- do.call("dgp_eta", specs[[j - 1]])
    eta$u <- eta$eta0 + rnorm(nrow(eta))
    d[[i <- paste("cat", j, sep = "")]] <- eta
    maps[[i]] <- attr(d[[i]], "maps")
  }
  d <- as.data.frame(d)

  u <- d[, grepl(".u", names(d), fixed = TRUE), drop = FALSE]
  u <- cbind("cat0.u" = 0, u)
  d$cat <- factor(apply(u, 1, function(x) { which(x == max(x)) }))

  for(j in 1:(nlevels - 1))
    names(d) <- gsub(paste("cat", j, ".x", sep = ""), "x", names(d))
  d <- d[, unique(names(d)), drop = FALSE]
  d <- d[, !grepl(".u", names(d))]

  attr(d, "maps") <- maps
  d
}


##########################
## (4) Simulation tests ##
##########################
if(FALSE) {
  ## Example 1:
  ## 3 category 2 level model.
  ## Create the data set.
  d <- dgp_multinomial(nlevels = 3, nobs = c(2000, 400),
    specs = list(
      list(
        type = list(c("sinus", "quadratic"), "double"),
        scale = list(c(1, 1), 1),
        const = 0
      ),
      list(
        type = list(c("sinus", "linear"), "double"),
        scale = list(c(-1, 1), -1),
        const = 0
      )
  ))

  ## The model formula, needs category specific id's!
  f <- list(
    cat ~ -1 + sx(x11) + sx(x12) + sx(cat2.id2, bs = "re"), cat2.id2 ~ sx(x21),
    cat ~ -1 + sx(x11) + sx(x12) + sx(cat3.id2, bs = "re"), cat3.id2 ~ sx(x21)
  )

  ## Start BayesX sampler.
  ## Set grid = NA, so fitted effects will be returned with original observations.
  b0 <- bayesr(f, data = d, family = multinomial, reference = 1,
    n.iter = 12000, burnin = 2000, thin = 10, verbose = TRUE, grid = NA,
    engine = "BayesX")

  ## Plot. estimated effects.
  plot(b0)

  ## Extract fitted values for e.g. term x11.
  fb0 <- fitted(b0, term = "x11")

  ## Plot estimated vs. true function
  ## Category 2
  plot2d(fb0[["2"]][[1]][[1]])
  with(d, lines(x11[order(x11)], cat2.f11[order(x11)], col = "red"))

  ## Category 3
  plot2d(fb0[["3"]][[1]][[1]])
  with(d, lines(x11[order(x11)], cat3.f11[order(x11)], col = "red"))

  ## Same with predict()
  nd <- d[, c("x11", "cat2.f11", "cat3.f11")]
  nd$p2 <- predict(b0, nd, model = c("2", "h1"), term = "x11")
  nd$p3 <- predict(b0, nd, model = c("3", "h1"), term = "x11")

  plot(p2[order(x11)] ~ x11[order(x11)], type = "l", data = nd)
  lines(cat2.f11[order(x11)] ~ x11[order(x11)], col = "red", data = nd)


  ## Example 2:
  ## Complicated spatial example, does not work for cat2 h2 spatial?!
  d <- dgp_multinomial(nlevels = 3, nobs = c(2000, 400),
    specs = list(
      list(
        type = list(c("sinus", "quadratic"), c("double", "spatial")),
        scale = list(c(1, 1), 1),
        const = 0
      ),
      list(
        type = list(c("sinus", "linear"), c("double", "spatial")),
        scale = list(c(-1, 1), -1),
        const = 0
      )
  ))

  ## Get the map and neighborhood matrix
  maps <- attr(d, "maps")
  nmat <- maps$cat2$map2$nmat
  pmap <- maps$cat2$map2$map

  ## Need to create new id's for the spatial effects, too.
  d$cat2.id20 <- d$cat2.id2
  d$cat3.id20 <- d$cat3.id2

  ## Model formula.
  f <- list(
    cat ~ -1 + sx(x11) + sx(x12) + sx(cat2.id2, bs = "re"),
      cat2.id2 ~ sx(x21) + sx(cat2.id20, bs = "mrf", map = nmat),
    cat ~ -1 + sx(x11) + sx(x12) + sx(cat3.id2, bs = "re"),
      cat3.id2 ~ sx(x21) + sx(cat2.id20, bs = "mrf", map = nmat)
  )

  ## Start BayesX sampler
  b1 <- bayesr(f, data = d, family = multinomial, reference = 1,
    n.iter = 12000, burnin = 2000, thin = 10, verbose = TRUE,
    engine = "BayesX")

  ## Plot effects with corresponding pixel map.
  plot(b1, map = pmap, legend = FALSE, scale = 0)
}

