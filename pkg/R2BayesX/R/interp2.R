interp2 <- function(x, y, z, xo = NULL, yo = NULL, k = 40, grid = 30, extrap = TRUE)
{
  if(is.null(xo))
    xo <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length = grid)
  if(is.null(yo))
    yo <- seq(min(y, na.rm = TRUE), max(y, na.rm = TRUE), length = grid)

  grid <- length(xo)

  x <- as.numeric(x); y <- as.numeric(y); z <- as.numeric(z)
  xo <- as.numeric(xo); yo <- as.numeric(yo)

  xr <- range(x, na.rm = TRUE)
  yr <- range(y, na.rm = TRUE)

  x <- (x - xr[1]) / diff(xr)
  y <- (y - yr[1]) / diff(yr)

  b <- mgcv::gam(z ~ s(x, y, k = k))

  x2 <- (xo - xr[1]) / diff(xr)
  y2 <- (yo - yr[1]) / diff(yr)

  nd <- data.frame("x" = rep(x2, grid), "y" = rep(y2, rep(grid, grid)))

  fit <- as.vector(predict(b, newdata = nd))

  if(!extrap) {
    require("sp")
    pid <- chull(X <- cbind(x, y))
    pol <- X[c(pid, pid[1]), ]
    pip <- point.in.polygon(nd$x, nd$y, pol[, 1], pol[, 2])
    fit[!pip] <- NA
  }

  return(matrix(fit, grid, grid))
}

