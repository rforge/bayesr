get.kde.range <- function(x, ngrid = 100L, probs = c(0.05, 0.95)) 
{
  kde <- stats::density(x, n = length(x), from = min(x), to = max(x))
  f <- stats::splinefun(x = kde$x, y = kde$y)
  xgrid <- matrix(seq(min(x), max(x), length = ngrid), ncol = 1L)
  foo <- function(f, lower, x) integrate(f, lower = lower, upper = x)$value
  xgridint <- apply(xgrid, 1L, foo, f = f, lower = min(x))
  q01 <- xgrid[xgridint >= probs[1L]]
  q01 <- q01[1L]
  q99 <- xgrid[xgridint <= probs[2L]]
  q99 <- q99[length(q99)]

  return(c(q01, q99))
}


get.kde.range2 <- function(x, ngrid = 100L, probs = c(0.05, 0.95)) 
{
  ex <- stats::ecdf(x)
  F <- ex(x)
  q01 <- x[F >= probs[1L]]
  q01 <- q01[1L]
  q02 <- x[F >= probs[2L]]
  q02 <- q02[1L]

  return(c(q01, q02))
}

