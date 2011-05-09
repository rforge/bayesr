plotmap <- function(map, x = NULL, id = NULL, c.select = NULL, legend = TRUE, 
  swap = FALSE, range = NULL, names = FALSE, values = FALSE, col = NULL, ncol = 100, 
  breaks = NULL, cex.legend = 1, cex.names = 1, cex.values = cex.names, digits = 2L,
  mar.min = 2, add = FALSE, ...)
{
  if(missing(map))
    stop("map object is missing!")
  if(!is.list(map))
    stop("argument map must be a list() of matrix polygons!")
  args <- list(...)
  map.limits <- find.limits(map, mar.min, ...)
  if(is.null(args$asp))
    args$asp <- attr(map, "asp")
  n <- length(map)
  if(inherits(x, "geo.bayesx"))
    names(map) <- as.character(1L:n)
  if(is.null(x))
    legend <- FALSE
  if(!any(is.na(poly.names <- f2int(names(map))))) {
      poly.names <- sort(poly.names)
      poly.names <- as.character(poly.names)
  } else poly.names <- sort(names(map))
  map <- map[poly.names]
  poly.names <- names(map)
  surrounding <- attr(map, "surrounding")
  inner.id <- which(sapply(surrounding, length) > 0L)
  if(length(inner.id)) {
    poly.names <- c(poly.names[- inner.id], poly.names[inner.id])
    map <- c(map[- inner.id], map[inner.id])
  }
  if(!is.null(args$ylim))
    map.limits$ylim <- args$ylim
  if(!is.null(args$xlim))
    map.limits$xlim <- args$xlim
  if(is.null(args$symmetric))
    symmetric <- TRUE
  else
    symmetric <- args$symmetric
  if(!is.null(x)) {
    if(is.null(col))
      col <- colorspace::diverge_hcl
    x <- compute.x.id(x, id, c.select, range, symmetric)
    colors <- make_pal(col = col, ncol = ncol, data = x$x, 
      range = range, breaks = breaks, swap = swap, 
      symmetric = symmetric)$map(x$x)
  } else {
    if(is.null(col))
      colors <- rep(NA, length.out = n)
    else {
      if(is.function(col))
        colors <- col(ncol)
      else colors <- col
      colors <- rep(colors, length.out = n)
    }
  }
  if(!is.null(map.limits$mar) && is.null(args$asp))
    par(mar = map.limits$mar)
  args$x <- map.limits$x
  args$y <- map.limits$y
  if(is.null(args$type))
    args$type <- "n"
  if(is.null(args$axes))
    args$axes <- FALSE
  if(is.null(args$xlab))
    args$xlab <- ""
  if(is.null(args$ylab))
    args$ylab <- ""
  if(!add)
    do.call(graphics::plot.default, delete.args(graphics::plot.default, args)) 
  args$ylab <- args$xlab <- args$main <- ""
  args$type <- NULL
  args$axes <- NULL
  lwd.p <- rep(args$lwd, length.out = n)
  if(is.null(lwd.p))
    lwd.p <- rep(0.2, length.out = n)
  lty.p <- rep(args$lty, length.out = n)
  if(is.null(lty.p))
    lty.p <- rep(1, length.out = n)
  border.p <- rep(args$border, length.out = n)
  if(is.null(border.p))
    border.p <- rep("black", length.out = n)
  density.p <- rep(args$density, length.out = n)
  angle.p <- rep(args$angle, length.out = n)
  if(is.null(angle.p))
    angle.p <- rep(45, length.out = n)
  i <- 1L
  for(poly in poly.names) {
    args$x <- map[[poly]][,1L]
    args$y <- map[[poly]][,2L]
    args$border <- border.p[i]
    args$angle <- angle.p[i]
    args$lwd <- lwd.p[i]
    args$lty <- lty.p[i]
    if(!is.null(density.p))
      args$density <- density.p[i]
    if(!is.null(x)){ 
      if(!is.na(k <- pmatch(poly, x$id))) {
        args$col <- colors[k]
        args$density <- NULL
      } else args$col <- NULL
    } else args$col <- colors[i]
    do.call(graphics::polygon, 
    delete.args(graphics::polygon, args, 
    c("lwd", "cex", "lty")))
    if(names && !values) {
      args$polygon <- map[[poly]]
      args$poly.name <- poly
      args$counter <- i
      args$cex <- cex.names
      do.call(centroidtext, delete.args(centroidtext, args, "font"))
    }
    if(values && !names) {
      args$polygon <- map[[poly]]
      args$poly.name <- as.character(round(x$x[k], digits = digits))
      args$counter <- k
      args$cex <- cex.values
      do.call(centroidtext, delete.args(centroidtext, args, "font"))
    }
    i <- i + 1L
  }
  if(legend) {
    if(is.null(args$pos))
      args$pos <- "bottomleft"
    if(is.null(args$width))
      args$width <- 0.4
    if(is.null(args$height))
      args$height <- 0.12
    if(is.null(args$distance.labels))
     args$distance.labels <- 1.5
    if(is.null(args$length.ticks))
      args$length.ticks <- 2L
    args$xlim <- map.limits$xlim
    args$ylim <- map.limits$ylim
    args$color <- col
    args$ncol <- ncol
    args$x <- x$x
    args$breaks <- breaks
    args$swap <- swap
    args$plot <- FALSE
    args$digits <- digits
    args$cex.labels <- cex.legend
    if(is.null(args$xpd))
      args$xpd <- TRUE
    args$symmetric <- symmetric
    args$range <- range
    args$add <- TRUE
    do.call(colorlegend, delete.args(colorlegend, args, c("font")))
  }
  if(!is.null(args$xlab))
    mtext(args$xlab, side = 1L)
  if(!is.null(args$ylab))
    mtext(args$ylab, side = 2L)

  return(invisible(NULL))
}


compute.x.id <- function(x, id, c.select, range, symmetric)
{
  if(is.null(id) && (is.vector(x) || is.array(x))) {
    if(!is.null(names(x))) {
      id <- names(x)
      x <- as.vector(x)
    }
  }
  if(is.factor(id))
    id <- f2int(id)
  if(is.array(x) && length(dim(x)) < 2L)
    x <- as.vector(x)
  if(is.vector(x) && is.vector(id)) {
    if(length(x) != length(id))
      stop("arguments x and id are differing!")
  } else {
    x <- unclass(x)
    if(is.list(x)) 
      nx <- names(x)
    if(is.matrix(x)) {
      x <- as.list(as.data.frame(x))
      nx <- names(x)  
      if(all(nx %in% paste("V", 1L:length(nx), sep = ""))) {
        nx[1L:2L] <- c("id", "x")
        c.select <- "x"
      }
    }
    if(is.data.frame(x)) {
      x <- as.list(x)
      nx <- names(x)
    }
    if(is.null(id))
      id <- x[[1L]]
    else {
      if(is.character(id)) {
        if(is.na(id <- pmatch(id, nx)))
          stop("argument id is specified wrong!")
      } else {
        if(id > length(nx))
          stop("argument id is specified wrong!")
      }
      id <- x[[id]]
    }
    if(is.null(c.select)) {
      take <- c("mean", "Mean", "MEAN", "estimate", 
        "Estimate", "ESTIMATE", "mean", "pmode", "pmean_tot")
      for(k in take)
        if(!is.na(pmatch(k, nx)))
          x <- x[[k]]
    } else {
      if(is.character(c.select)) {
        k <- pmatch(c.select, nx)
      if(is.na(k))
        stop("argument c.select is specified wrong!")
      x <- x[[k]]
      } else {
        if(c.select > length(nx))
          stop("argument c.select is specified wrong!")
        x <- x[[c.select]]
      }
    }
  }
  if(symmetric) {
    if(is.null(range)) {
      if(min(x) < 0)
        m <- (-1)
      else
        m <- 1
      if(abs(min(x)) > abs(max(x)))
        x <- c(x, abs(min(x)))
      if(abs(max(x)) > abs(min(x)))
        x <- c(x, m * abs(max(x)))
      id <- c(as.character(id), "added")
    } else {
      if(max(range) > max(x)) {
        x <- c(x, max(range))
        id <- c(as.character(id), "added")
      } else x[x > max(range)] <- max(range)
      if(min(range) < min(x)) {
        x <- c(x, min(range))
        id <- c(as.character(id), "added")
      } else x[x < min(range)] <- min(range)
    }
  }

  return(list(id = as.character(id), x = x))
}


find.limits <- function(map, mar.min = 2, ...)
{
  if(!is.list(map))
    stop("argument map must be a list() of matrix polygons!")
  n <- length(map)
  myrange <- function(x, c.select = 1L, ...) {
    return(na.omit(x[,c.select], ...))
  }
  xlim <- range(unlist(lapply(map, myrange, c.select = 1L, ...)))
  ylim <- range(unlist(lapply(map, myrange, c.select = 2L, ...)))
  mar <- NULL
  if(!is.null(height2width <- attr(map, "height2width"))) {
    height2width <- height2width * 1.1
    if(!is.null(mar.min)) {
      if(height2width > 1) {
        side <- 17.5*(1-1/height2width)+mar.min/height2width
        mar <- c(mar.min, side, mar.min, side)
      }
      else {
        top <- 17.5*(1-height2width)+mar.min*height2width
        mar <- c(top, mar.min, top, mar.min)
      }
    }
  }

  return(list(ylim = ylim, xlim = xlim, mar = mar))
}
