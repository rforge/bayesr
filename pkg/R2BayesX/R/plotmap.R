plotmap <- function(map, x = NULL, id = NULL, c.select = NULL, legend = TRUE,
  swap = FALSE, range = NULL, names = FALSE, values = FALSE, col = NULL,
  ncol = 100, breaks = NULL, cex.legend = 1, cex.names = 1, cex.values = cex.names,
  digits = 2L, mar.min = 2, add = FALSE, interp = FALSE, linear = FALSE, extrap = FALSE,
  outside = FALSE, grid = 100, p.pch = 16, p.cex = 1, ...)
{
  if(missing(map))
    stop("map object is missing!")
  if(inherits(map, "SpatialPolygons"))
    map <- sp2bnd(map)
  if(!is.list(map))
    stop("argument map must be a list() of matrix polygons!")
  args <- list(...)
  map.limits <- find.limits(map, mar.min, ...)
  if(is.null(args$asp))
    args$asp <- attr(map, "asp")
  n <- length(map)
  if(is.null(x))
    legend <- FALSE
  poly.names.orig <- names(map)
  if(!any(is.na(poly.names <- x2int(names(map))))) {
    poly.names <- sort(poly.names)
    poly.names <- as.character(poly.names)
  } else {
    poly.names <- sort(names(map))
  }
  if(length(upn <- unique(poly.names)) < length(poly.names)) {
    nn <- NULL
    for(i in upn) {
      j <- poly.names == i
      poly.names[j] <- paste(poly.names[j],
        if(sum(j) > 1) paste(".", 1:sum(j), sep = "") else NULL,
        sep = "")
    }
    names(map) <- poly.names
  }
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
    if(is.null(col)) {
      col <- colorspace::diverge_hcl
      # col <- colorspace::diverge_hcl(ncol, h = c(130, 10), c = 250,
      #  l = c(30, 90), power = 1.5, gamma = 2.4, fixup = TRUE)
    }
    x <- compute.x.id(x, id, c.select, range, symmetric)
    map_fun <- make_pal(col = col, ncol = ncol, data = x$x, 
      range = range, breaks = breaks, swap = swap, 
      symmetric = symmetric)$map
    colors <- map_fun(x$x)
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
  if(legend && !is.null(args$pos) && args$pos[1L] == "right") {
    par.orig <- par(c("mar", "las", "mfrow"))
    mar.orig <- mar <- par.orig$mar
    mar[4L] <- 0
    on.exit(par(par.orig))
    par(mar = mar)
    w <- (3 + mar[2L]) * par("csi") * 2.54
    layout(matrix(c(1, 2), nrow = 1), widths = c(1, lcm(w)))
  }
  if(!is.null(map.limits$mar) && is.null(args$asp) && !add)
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
  if(interp & !is.null(x)) {
    stopifnot(require("maptools"))

    cdata <- data.frame(R2BayesX:::centroids(map), "id" = names(map))
    cdata <- merge(cdata, data.frame("z" = x$x, "id" = x$id), by = "id")

    ico <- with(cdata, akima::interp(x = xco, y = yco, z = z,
      xo = seq(map.limits$x[1], map.limits$x[2], length = grid),
      yo = seq(map.limits$y[1], map.limits$y[2], length = grid),
      duplicate = "strip", linear = linear, extrap = extrap))
    
    yco <- rep(ico$y, each = length(ico$x))
    xco <- rep(ico$x, length(ico$y))

    cvals <- as.numeric(ico$z)
    cvals[cvals < min(cdata$z)] <- min(cdata$z)
    cvals[cvals > max(cdata$z)] <- max(cdata$z)
    icolors <- map_fun(cvals)

    if(!outside) {
      gpclibPermit()
      class(map) <- "bnd"
      mapsp <- bnd2sp(map)
      ob <- unionSpatialPolygons(mapsp, rep(1L, length = length(mapsp)), avoidGEOS  = TRUE)

      nob <- length(slot(slot(ob, "polygons")[[1]], "Polygons"))
      pip <- NULL
      for(j in 1:nob) {
        oco <- slot(slot(slot(ob, "polygons")[[1]], "Polygons")[[j]], "coords")
        pip <- cbind(pip, point.in.polygon(xco, yco, oco[, 1L], oco[, 2L], mode.checked = FALSE) < 1L)
      }
      pip <- apply(pip, 1, function(x) all(x))
    
      icolors[pip] <- NA
    }

    points(SpatialPoints(cbind(xco, yco)), col = icolors, pch = p.pch, cex = p.cex)
    colors <- rep(NA, length = length(colors))
  }
  args$ylab <- args$xlab <- args$main <- ""
  args$type <- NULL
  args$axes <- NULL
  lwd.p <- rep(args$lwd, length.out = n)
  if(is.null(lwd.p))
    lwd.p <- rep(1, length.out = n)
  lty.p <- rep(args$lty, length.out = n)
  if(is.null(lty.p))
    lty.p <- rep(1, length.out = n)
  border.p <- rep(args$border, length.out = n)
  if(is.null(border.p))
    border.p <- rep("black", length.out = n)
  density.p <- rep(args$density, length.out = n)
  angle.p <- rep(args$angle, length.out = n)
  if(is.null(angle.p))
    angle.p <- rep(90, length.out = n)
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
      } else {
        args$col <- NULL
        if(is.null(args$density))
          args$density <- 20L
      }
    } else args$col <- colors[i]
    do.call(graphics::polygon, 
    delete.args(graphics::polygon, args, 
    c("lwd", "cex", "lty")))
    if(names && !values) {
      args$polygon <- map[[poly]]
      args$poly.name <- poly.names.orig[i]
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
      args$pos <- "topleft"
    if(args$pos[1L] == "right") {
      args$full <- TRUE
      args$side.legend <- 2L
      args$side.ticks <- 2L
      mar <- mar.orig
      mar[2L] <- 0.5
      mar[4L] <- 3.1
      par(mar = mar, xaxs = "i", yaxs = "i")
      args$plot <- TRUE
      args$add <- FALSE
    } else {
      args$plot <- FALSE
      if(is.null(args$xpd))
        args$xpd <- TRUE
      args$add <- TRUE
    }
#    if(is.null(args$width))
#      args$width <- 0.6
#    if(is.null(args$height))
#      args$height <- 0.2
    if(is.null(args$distance.labels))
      args$distance.labels <- 1
    if(is.null(args$length.ticks))
      args$length.ticks <- 2L
    args$xlim <- map.limits$xlim
    args$ylim <- map.limits$ylim
    args$color <- col
    args$ncol <- ncol
    args$x <- x$x
    args$breaks <- breaks
    args$swap <- swap
    args$digits <- digits
    args$cex.labels <- cex.legend
    args$symmetric <- symmetric
    if(is.null(range)) {
      range <- range(args$x)
      if(diff(range) == 0)
        range <- unique(range) + c(-1, 1)
    }
    args$range <- range
    if(is.null(args$lrange))
      args$lrange <- args$range
    do.call(colorlegend, delete.args(colorlegend, args, c("font")))
  }
  if(!is.null(args$xlab))
    mtext(args$xlab, side = 1L)
  if(!is.null(args$ylab))
    mtext(args$ylab, side = 2L)

  return(invisible(NULL))
}

