GAMart <- function(n = 500, sd = 0.6, seed = TRUE)
{
  if(seed) set.seed(111)

  ## help scale function
  hs <- function(x, min = 0.1, max = 0.6) {
    x <- if(length(unique(x)) > 1) {
      (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE)) * (max - min) + min
    } else x
    x
  }
     
  ## (1) regressors
  ## spatial
  n2 <- ceiling(sqrt(n) / 4)
  d <- expand.grid("long" = seq(0, 1, length = n2), "lat" = seq(0, 1, length = n2))
  d <- rep(d, ceiling(n / (n2 * n2)))
  d <- data.frame("long" = unlist(d[grepl("long", names(d))]),
    "lat" = unlist(d[grepl("lat", names(d))]))
  d <- d[1:n, ]
  
  ## other covariates
  d$x1 <- runif(n, 0, 1)
  d$x2 <- runif(n, 0, 1)
  d$x3 <- runif(n, 0, 1)

  ## (4) functions
  ## linear
  f1 <- function(x) {
    hs(-1.2 * x, -1, 1)
  }

  ## doublemode
  f2 <- function(x) {
    hs(1.3 * (120 * x * exp(-10 * x) + 2.75 * x^2), -1, 1)
  }

  ## quadratic 
  f3 <- function(x) {
    hs(3.5 * (x - mean(x))^2, -1, 1)
  }

  ## spatial
  f4 <- function(long, lat) {
    hs(sin(hs(long, -3, 3)) * cos(hs(lat, -3, 3)), -1, 1)
  }
  
  ## response
  d$f1 <- f1(d$x1); d$f2 <- f2(d$x2); d$f3 <- f3(d$x3); d$f4 <- f4(d$long, d$lat);
  d$eta <- with(d, hs(f1 + f2 + f3 + f4, -1, 1))
  d$num <- with(d, eta + rnorm(n, sd = sd))
  d$bin <- 1 * (hs(d$eta, -1, 1) + rnorm(n) > 0)
  d$cat <- cut(d$num, quantile(d$num), include.lowest = TRUE)

  d
}

