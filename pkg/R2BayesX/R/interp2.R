interp2 <- function(x, y, z, xo = NULL, yo = NULL, grid = 30,
  type = c("akima", "mba", "mgcv"), linear = TRUE, extrap = FALSE, k = 40)
{
  type <- tolower(type)
  type <- match.arg(type)

  if(is.null(xo))
    xo <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length = grid)
  if(is.null(yo))
    yo <- seq(min(y, na.rm = TRUE), max(y, na.rm = TRUE), length = grid)

  grid <- length(xo)
  x <- as.numeric(x); y <- as.numeric(y); z <- as.numeric(z)

  if(type == "mgcv") {
    xo <- as.numeric(xo); yo <- as.numeric(yo)
    xr <- range(x, na.rm = TRUE)
    yr <- range(y, na.rm = TRUE)
    x <- (x - xr[1]) / diff(xr)
    y <- (y - yr[1]) / diff(yr)

    if(k > length(z))
      k <- ceiling(0.8 * length(z))

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
  }
  if(type == "mba") {
    stopifnot(require("MBA"))
    fit <- MBA::mba.surf(data.frame("x" = x, "y" = y, "z" = z), grid, grid)$xyz.est$z
  }
  if(type == "akima") {
    if(isTRUE(getOption("use.akima"))) {
      stopifnot(require("akima"))
    } else {
      if(require("akima")) {
        cat("NOTE: Package 'akima' has an ACM license that restricts applications to non-commercial usage.\n")
      } else {
        stop(paste("plot3() can only be used if the 'akima' package is installed. ",
          "Note that 'akima' has an ACM license that restricts applications to ",
          "non-commercial usage.", sep = ""))
      }
    }

    fit <- try(akima::interp(x, y, z, xo = xo, yo = yo, 
      duplicate = "strip", linear = linear, extrap = extrap)$z, silent = TRUE)
    if(inherits(fit, "try-error")) {
      fit <- try(akima::interp(jitter(x), jitter(y), z, xo = xo, yo = yo, 
        duplicate = "strip", linear = linear, extrap = extrap)$z, silent = TRUE)
    }
  }

  return(matrix(fit, grid, grid))
}

