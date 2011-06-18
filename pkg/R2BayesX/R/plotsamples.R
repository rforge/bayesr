plotsamples <-
function(x, selected = "NA", acf = FALSE, var = FALSE, max = FALSE, ...)
{
  if(is.null(x)) {
    warning("there is nothing to plot!")
    return(invisible(NULL))
  }
  args <- list(...)
  if(is.null(args$axes))
    axes <- TRUE
  else
    axes <- args$axes
  xlab <- args$xlab
  ylab <- args$ylab
  main <- args$main
  if(!is.matrix(x))
    x <- matrix(x, length(x), 1L)
  nr <- ncol(x)
  op <- par(no.readonly = TRUE)
  args$xlab <- NA
  args$ylab <- NA
  if(is.null(args$main))
    args$main <- NA
  if(!acf) {
    args$type <- "l"
    args$x <- 1L:nrow(x)
  } else args$verbose <- FALSE
  if(nr > 1L && !max)
    setmfrow(nr)
  if(var) {
    outer <- FALSE
    line <- 2L
  } else {
    if(!max) {
      if(nr < 2L)
        par(mar = c(2, 2, 2, 2))
      else
        par(mar = c(2, 2, 3, 2))
    }
    outer <- TRUE
    line <- 0L
  }
  if(max)
    maxs <- NULL
  for(k in 1L:nr) {
    if(!acf) {
      args$y <- x[,k]
      args$axes <- FALSE
      if(max)
        maxs <- cbind(maxs, x[,k])
      else
        do.call(graphics::plot.default, args)
    } else {
      args$x <- stats::ts(data = x[,k])
      args$plot <- FALSE
      acfx <- do.call(stats::acf, args)
      lx <- as.vector(acfx$lag)
      ax <- as.vector(acfx$acf)
      n <- length(lx)
      lx <- lx[2L:n]
      ax <- ax[2L:n]
      acfx$lag <- array(lx, dim = c(n - 1L, 1L, 1L))
      acfx$acf <- array(ax, dim = c(n - 1L, 1L, 1L))
      if(max) {
        maxs <- cbind(maxs, acfx$acf)
      } else {
        ylim <- NULL
        if(is.null(args$ylim) && !var)
          ylim <- c(-0.2, 1)
        else
          ylim = args$ylim
        if(is.null(args$main))
          args$main <- NA
        stats:::plot.acf(acfx, main = args$main, axes = FALSE, ylim = ylim)
      }
    }
    if(!max) {
      if(axes) {
        box()
        at <- axis(1L, tick = FALSE, labels = NA)
        at[1L] <- 1L
        axis(1L, at = at, labels = at)
        axis(2L)
      }
      if(nr > 1L)
        mtext(paste("Coefficient", k), side = 3L, line = 0.5, cex = 0.8)
      if(!is.null(xlab))
        mtext(xlab, side = 1L, line = line, outer = outer, font = 1L)
      if(!is.null(ylab))
        mtext(ylab, side = 2L,line = line, outer = outer, font = 1L)
      if(!is.null(main))
        mtext(main, side = 3L, line = line, outer = outer, font = 2L, cex = 1)
      if(is.null(main)) {
        if(var)
          ptxt <- "Variance"
        else
          ptxt <- "Coeffiecient(s)"
        if(acf)
          ptxt <- paste(ptxt, "autocorrelation")
        else
          ptxt <- paste(ptxt, "sampling path(s)")
        ptxt <- paste(ptxt, "of term", selected)
        mtext(ptxt, side = 3L, line = line, outer = outer, font = 2L, cex = 1)
      }
    }
  } 
  if(max) {
    maxs <- apply(maxs, 1L, max)
    if(acf) {
      n <- length(maxs)
      acfx$lag <- array(lx, dim = c(n, 1L, 1L))
      acfx$acf <- array(maxs, dim = c(n, 1L, 1L))
      ylim <- NULL
      if(is.null(args$ylim) && !var)
        ylim <- c(-0.2, 1)
      else
        ylim = args$ylim
      if(is.null(args$main) || is.na(args$main)) {
        acfx$main <- "Maximum autocorrelation of parameters"
        if(selected != "NA" && !is.null(selected))
          acfx$main <- paste(acfx$main, "of term", selected)
      }
      else
        acfx$main <- args$main
      stats:::plot.acf(acfx, main = acfx$main, axes = FALSE, ylim = ylim)
    } else {
      args$y <- maxs
      if(is.null(args$main) || is.na(args$main)) {
        args$main <- "Maximum of samples of parameters"
        if(selected != "NA" && !is.null(selected))
          args$main <- paste(args$main, "of term", selected) 
      }
      if(is.null(args$xlab) || is.na(args$xlab))
        args$xlab <- "Iteration"
      if(is.null(args$ylab) || is.na(args$xlab))
        args$ylab <- ""
      if(is.null(args$type))
        args$type <- "l"
      args$acf <- NULL
      args$lag.max <- NULL
      do.call(graphics::plot.default, args)
    }
    if(axes) {
      box()
      at <- axis(1L, tick = FALSE, labels = NA)
      at[1L] <- 1L
      axis(1L, at = at, labels = at)
      axis(2L)
    }
  }
  if(nr > 1L)
    par(op)

  return(invisible(NULL))
}

