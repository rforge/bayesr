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
scale3 <- function(x, lower = -1.5, upper = 1.5)
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
        fun[[paste("sp", j, sep = "")]] <- scale3(e$f) * if(is.na(scale[[j]][i])) 1 else scale[[j]][i]
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
          scale3(FUN(x)) * if(is.na(scale[[j]][i])) 1 else scale[[j]][i]
        } else FUN(const = const)
        if(mt != "const") xd[[xn]] <- x
      }
    }
    eta <- 0
    for(i in seq_along(fun)) {
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
  sd <- scale3(sigma$eta0, range.sigma[1], range.sigma[2])
  d$y <- rnorm(nrow(mu), mean = mu$eta0, sd = sd)

  d
}

if(FALSE) {
  d <- dgp_gaussian()
  b <- bayesx2(y ~ sx(mu.x11) + sx(mu.x12), ~ sx(sigma.x11), data = d)
  d$p <- predict(b, model = "mu", term = c("x11", "x12"))
}


## Gamma.
dgp_gamma <- function(mu = NULL, sigma = NULL, ...)
{

}


## Beta.
dgp_beta <- function(mu = NULL, sigma = NULL, ...)
{

}


## Multivariate normal.
dgp_mvn <- function(n = 500, mu1 = NULL, mu2 = NULL, sigma1 = NULL, sigma2 = NULL, rho = NULL,
  range.sigma1 = c(0.3, 1.2), range.sigma2 = c(0.3, 1.2), range.rho = c(0.3, 1.2), ...)
{
  require("mvtnorm")

  if(is.null(mu1)) {
    mu1 <- list(nobs = n, const = 0.5,
      type = list(c("unimodal", "quadratic", "const")))
  }
  if(is.null(mu2)) {
    mu2 <- list(nobs = n, const = 0.1,
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
    sm <- cbind(c(s1, p * s1 * s2), c(p * s1 * s2, s2))
    y <- rbind(y, rmvnorm(1, m, sm, method = "svd"))
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
    sigma1 ~ sx(sigma1.x11),
    sigma2 ~ sx(sigma2.x11),
    rho ~ sx(rho.x11)
  )

  b <- bayesr(f, family = mvn, data = d, engine = "BayesX", verbose = TRUE)
}


## Multinomial.
dgp_multinom <- function(nlevels = 4,
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
  d <- dgp_multinom(nlevels = 3, nobs = c(2000, 400),
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
  d <- dgp_multinom(nlevels = 3, nobs = c(2000, 400),
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

