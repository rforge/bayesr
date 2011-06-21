plotsamples <-
function(x, selected = "NA", acf = FALSE, var = FALSE, all.acf = FALSE, ...)
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
  if(!acf && !all.acf)
    args$lag.max <- NULL
  xlab <- args$xlab
  ylab <- args$ylab
  main <- args$main
  if(!is.matrix(x))
    x <- matrix(x, length(x), 1L)
  nr <- ncol(x)
  op <- par(no.readonly = TRUE)
  if(is.null(args$xlab)) {
    if(acf || all.acf)
      args$xlab <- "Lag"
    else
      args$xlab <- "Iteration"
  }
  if(is.null(args$ylab)) {
    if(acf || all.acf)
      args$ylab <- "ACF"
    else
      args$ylab <- "Sample"
  }
  if(is.null(args$main))
    args$main <- NA
  if(!acf && !all.acf) {
    args$type <- "l"
    args$x <- 1L:nrow(x)
  } else args$verbose <- FALSE
  if(nr > 1L && !all.acf)
    setmfrow(nr)
  line <- 2L
  sa <- TRUE
  if(nr < 2L || all.acf)
    outer <- FALSE
  else {
    if(is.null(args$oma) && nr > 2L && all(par()$oma == c(0, 0, 0, 0)))
      par(oma = c(5.1, 4.1, 5.1, 2.1))
    if(is.null(args$mar) && nr > 2L && all(par()$mar == c(5.1, 4.1, 4.1, 2.1)))
      par(mar = c(4.8, 4.1, 2.1, 1.5))
    outer <- TRUE
  }
  if(all.acf)
    maxs <- NULL
  for(k in 1L:nr) {
    if(!acf && !all.acf) {
      args$y <- x[,k]
      args$axes <- FALSE
      do.call(graphics::plot.default, args)
    } else {
      args$x <- stats::ts(data = x[,k])
      args$plot <- FALSE
      if(is.null(args$lag.max))
        args$lag.max <- length(x[,k])
      acfx <- do.call(stats::acf, args)
      lx <- as.vector(acfx$lag)
      ax <- as.vector(acfx$acf)
      n <- length(lx)
      lx <- lx[2L:n]
      ax <- ax[2L:n]
      acfx$lag <- array(lx, dim = c(n - 1L, 1L, 1L))
      acfx$acf <- array(ax, dim = c(n - 1L, 1L, 1L))
      if(all.acf) {
        maxs <- cbind(maxs, acfx$acf)
      } else {
        ylim <- NULL
        if(is.null(args$ylim))
          ylim <- c(-0.2, 1)
        else
          ylim = args$ylim
        stats:::plot.acf(acfx, main = args$main, axes = FALSE, ylim = ylim, xlab = args$xlab,
          ylab = args$ylab)
      }
    }
    if(!all.acf) {
      if(axes) {
        box()
        at <- axis(1L, tick = FALSE, labels = NA)
        at[1L] <- 1L
        axis(1L, at = at, labels = at)
        axis(2L)
      }
      if(nr > 1L) {
        if(var)
          mtext(paste("Parameter", k), side = 3L, line = 0.5, cex = 0.8)
        else
          mtext(paste("Coefficient", k), side = 3L, line = 0.5, cex = 0.8)
      }
      if(!is.null(xlab) && !acf)
        mtext(xlab, side = 3L, line = line, outer = outer, font = 1L)
      if(!is.null(ylab) && !acf)
        mtext(ylab, side = 3L,line = line, outer = outer, font = 1L)
      if(!is.null(main) && !acf)
        mtext(main, side = 5L, line = line, outer = outer, font = 2L, cex = 1)
      if(is.null(main)) {
        if(var)
          ptxt <- "Variance"
        else {
          if(nr > 1L)
            ptxt <- "Coeffiecients"
          else
            ptxt <- "Coeffiecient"
        }
        if(acf)
          ptxt <- paste(ptxt, "autocorrelation")
        else {
          if(nr > 1L)
            ptxt <- paste(ptxt, "sampling paths")
          else
            ptxt <- paste(ptxt, "sampling path")
        }
        ptxt <- paste(ptxt, "of term", selected)
        mtext(ptxt, side = 3L, line = line, outer = outer, font = 2L, cex = 1)
      }
    }
  } 
  if(all.acf) {
    amax <- apply(maxs, 1L, max)
    amin <- apply(maxs, 1L, min)
    amax[abs(amin) > abs(amax)] <- amin[abs(amin) > abs(amax)]
    maxs <- amax
    n <- length(maxs)
    acfx$lag <- array(lx, dim = c(n, 1L, 1L))
    acfx$acf <- array(maxs, dim = c(n, 1L, 1L))
    ylim <- NULL
    if(is.null(args$ylim))
      ylim <- c(-0.2, 1)
    else
      ylim = args$ylim
    if(is.null(args$main) || is.na(args$main)) {
      acfx$main <- "Greatest upper/lower autocorrelation \nof sampled parameters"
      if(selected != "NA" && !is.null(selected))
        acfx$main <- paste(acfx$main, "of term", selected)
      if(all(par()$mar == c(5.1, 4.1, 4.1, 2.1)))
        par(mar = c(5.1, 4.1, 5.1, 2.1))
    } else acfx$main <- args$main
    stats:::plot.acf(acfx, main = acfx$main, axes = FALSE, ylim = ylim, xlab = args$xlab,
      ylab = args$ylab)
    if(axes) {
      box()
      at <- axis(1L, tick = FALSE, labels = NA)
      at[1L] <- 1L
      axis(1L, at = at, labels = at)
      axis(2L)
    }
  }
  if(nr > 1L && !all.acf)
    par(op)

  return(invisible(NULL))
}

