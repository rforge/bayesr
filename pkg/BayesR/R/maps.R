## Plot map using longitude and latitude coordinates
## using the maps database.
xymap <- function(x, y, z, color = sequential_hcl(99, h = 100), raw.color = FALSE, symmetric = FALSE,
  swap = TRUE, p.cex = 0.01, pch = 22, legend = TRUE, add = FALSE, domar = TRUE,
  layout = TRUE, mmar = c(0, 0, 0, 0), lmar = c(1.3, 3, 1.3, 2), interp = FALSE,
  grid = 8, linear = FALSE, extrap = FALSE, duplicate = "mean", xlim = NULL,
  ylim = NULL, boundary = TRUE, interior = TRUE, rivers = FALSE, mcol = NULL,
  contour.data = NULL, k = 30, akima = FALSE, data = NULL, subset = NULL, box = FALSE,
  ireturn = FALSE, sort = TRUE, ...)
{
  require("R2BayesX")
  require("maps")
  require("sp")

  if(!ireturn) {
    opar <- par(no.readonly = TRUE)
    on.exit(opar)
    if(domar & layout) {
      omar <- opar$mar
      par(mar = mmar)
    }
  }

  if(!is.null(data)) {
    subset <- deparse(substitute(subset), backtick = TRUE, width.cutoff = 500L)
    if(subset != "NULL") {
      subset <- eval(parse(text = subset), envir = data)
      data <- subset(data, subset)
    }
    x <- eval(parse(text = deparse(substitute(x), backtick = TRUE, width.cutoff = 500L)), envir = data)
    y <- eval(parse(text = deparse(substitute(y), backtick = TRUE, width.cutoff = 500L)), envir = data)
    z <- eval(parse(text = deparse(substitute(z), backtick = TRUE, width.cutoff = 500L)), envir = data)
  }

  data <- unique(data.frame("x" = as.numeric(x), "y" = as.numeric(y), "z" = as.numeric(z)))
  if(sort)
    data <- data[order(data$x), ]
  data <- na.omit(data)

  if(interp) {
    require("MBA")

    data <- na.omit(data)

    pp <- cbind(data$x, data$y)
    dx <- abs(diff(pp[, 1])); dy <- abs(diff(pp[, 2]))
    dx <- dx[dx != 0]; dy <- dy[dy != 0]
    res <- c(min(dx), min(dy)) / 2
    
    grid <- grid + 1

    px <- apply(pp, 1, function(x) {
      xs <- seq(x[1] - res[1], x[1] + res[1], length = grid)
      ys <- seq(x[2] - res[2], x[2] + res[2], length = grid)
      xs <- xs[-grid] + (xs[2] - xs[1]) / 2
      ys <- ys[-grid] + (ys[2] - ys[1]) / 2
      expand.grid("x" = xs, "y" = ys)
    })
    px <- do.call("rbind", px)

    data <- as.data.frame(mba.points(data, px, extend = TRUE, verbose = FALSE)$xyz.est)

    where <- map.where("world", data$x, data$y)
    data <- data[!is.na(where), ]
    if(ireturn)
      return(data)
  }
  
  colors <- colorlegend(x = data$z, plot = FALSE, color = color,
    swap = swap, symmetric = symmetric, ...)
  col <- colors$map(data$z)
  p.cex <- if(is.null(p.cex)) {
    if(interp) 0.1 else 1
  } else p.cex

  coordinates(data) <- c("x", "y")
  proj4string(data) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

  if(!add && legend && layout) {
    mar <- par()$mar
    par(mar = mar)
    w <- (4.5 + mar[4L]) * par("csi") * 3.7
  }
  if(!add && legend && layout)
    layout(matrix(c(1, 2), nrow = 1), widths = c(1, lcm(w)))

  pp <- coordinates(SpatialPoints(coordinates(data),
    proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")))
  dx <- abs(diff(pp[, 1])); dy <- abs(diff(pp[, 2]))
  dx <- dx[dx != 0]; dy <- dy[dy != 0]
  res <- c(min(dx), min(dy))

  if(!add) {
    plot(data, col = NA, bg = NA, xlim = xlim, ylim = ylim)
    if(box)
      box()
  }

  m <- map("world", add = TRUE, xlim = if(!add) NULL else xlim, ylim = if(!add) NULL else ylim,
    fill = if(is.null(mcol)) FALSE else TRUE, col = if(is.null(mcol)) gray(0.6) else mcol,
    boundary = boundary, interior = interior)

  if(rivers) {
    require("mapdata")
    map("rivers", add = TRUE, col = "lightblue")
  }

  rect(pp[, 1] - res[1] / 2, pp[, 2] - res[2] / 2, pp[, 1] + res[1] / 2, pp[, 2] + res[2] / 2,
    col = col, border = NA, lwd = 0)

  if(!is.null(mcol))
    points(data, col = col, bg = col, pch = pch, cex = p.cex)
  if(!is.null(contour.data)) {
    if(is.logical(contour.data) & contour.data)
      contour.data <- data.frame("x" = as.numeric(x), "y" = as.numeric(y), "z" = as.numeric(z))
    contour.data <- unique(contour.data)
    x <- contour.data[, 1]; y <- contour.data[, 2]; z <- contour.data[, 3]
    if(!akima) {
      cm <- bam(z ~ s(x, y, k = k))
      x <- seq(min(x), max(x), length = 100)
      y <- seq(min(y), max(y), length = 100)
      nd <- expand.grid(x = x, y = y)
      fit <- as.numeric(predict(cm, newdata = nd))
      ok <- map.where("world", nd$x, nd$y)
      fit[is.na(ok)] <- NA
      fit <- matrix(fit, 100, 100)
    } else {
      xo <- seq(min(x), max(x), length = 100)
      yo <- seq(min(y), max(y), length = 100)
      adat <- interp(x, y, z, xo = xo, yo = yo, duplicate = "mean")
      fit <- adat$z; x <- xo; y <- yo
    }
    contour(x, y, fit, add = TRUE)
  }
  if(legend) {
    if(layout) {
      par(mar = lmar)
      colorlegend(x = data$z, full = TRUE,
        side.legend = 2, side.ticks = 2, color = color,
        swap = swap, symmetric = symmetric, ...)
    } else {
      colorlegend(x = data$z, plot = FALSE, add = TRUE,
        color = color, swap = swap, symmetric = symmetric, ...)
    }
  }
  
  invisible(colors)
}


## Function to compute a polygon map from
## a grid overlayed on x- and y-coordinates.
pixelmap <- function(x, y, size = 0.1, width = NULL, data = NULL,
  all = TRUE, at.xy = FALSE, scale = TRUE, n = 20, ...)
{
  if(missing(x) & missing(y) & is.null(data)) {
    x <- expand.grid("x" = seq(0, 1, length = n), "y" = seq(0, 1, length = n))
    at.xy <- TRUE; size <- NA
  }
  if(missing(y) & (is.matrix(x) | is.data.frame(x) | is.list(x))) {
    x <- as.data.frame(x)
  } else {
    if(inherits(x, "formula")) {
      if(is.null(data))
        data <- environment(x)
      x <- model.frame(x, data = data)
      x <- x[, 2:1]
    } else {
      nx <- c(deparse(substitute(x), backtick = TRUE, width.cutoff = 500),
        deparse(substitute(y), backtick = TRUE, width.cutoff = 500))
      x <- data.frame(x, y)
      names(x) <- nx
    }
  }

  xr <- range(x[, 1], na.rm = TRUE)
  yr <- range(x[, 2], na.rm = TRUE)

  if(at.xy) {
    id <- rep(NA, nrow(x))
    xu <- unique(x)
    if(is.null(width)) {
      xd <- abs(diff(xu[, 1]))
      xd <- xd[xd > 0]
      xstep <- if(is.numeric(size)) {
        size * min(xd, na.rm = TRUE)
      } else {
        min(xd, na.rm = TRUE) / 2
      }
    } else xstep <- width / 2
    if(scale) {
      p <- xstep / abs(diff(xr))
      ystep <- abs(diff(yr)) * p
    } else ystep <- xstep
    map <- list()
    n <- nrow(xu)
    for(j in 1:n) {
      map[[j]] <- cbind(
        "x" = c(xu[j, 1] - xstep, xu[j, 1] + xstep, xu[j, 1] + xstep, xu[j, 1] - xstep, x[j, 1] - xstep),
        "y" = c(xu[j, 2] - ystep, xu[j, 2] - ystep, xu[j, 2] + ystep, xu[j, 2] + ystep, x[j, 2] - ystep)
      )
      id[x[, 1] >= xu[j, 1] - xstep & x[, 1] <= xu[j, 1] + xstep & x[, 2] >= xu[j, 2] - ystep & x[, 2] <= xu[j, 2] + ystep] <- j
    }
    names(map) <- as.character(1:n)
  } else {
    xstep <- if(is.null(width)) abs(diff(xr)) * size else width / 2
    if(scale) {
      p <- xstep / abs(diff(xr))
      ystep <- abs(diff(yr)) * p
    } else ystep <- xstep
    xstart <- xstart0 <- xr[1] - 0.5 * xstep
    ystart <- yr[1] - 0.5 * ystep
    map <- list(); k <- 1; id <- rep(NA, nrow(x))
    while(ystart < yr[2]) {
      xstart <- xstart0
      while(xstart < xr[2]) {
        map[[k]] <- cbind(
          "x" = c(xstart, xstart + xstep, xstart + xstep, xstart, xstart),
          "y" = c(ystart, ystart, ystart + ystep, ystart + ystep, ystart)
        )
        id[x[, 1] >= xstart & x[, 1] < xstart + xstep & x[, 2] >= ystart & x[, 2] < ystart + ystep] <- k
        xstart <- xstart + xstep
        k <- k + 1
      }
      ystart <- ystart + ystep
    }
    names(map) <- as.character(1:length(map))

    if(!all)
      map <- map[names(map) %in% as.character(id)]
  }

  class(map) <- c("bnd", "list")
  nmat <- neighbormatrix(map, ...)
  nn <- rowSums(nmat)
  nmat[nmat > 0] <- -1
  diag(nmat) <- nn
  id <- factor(as.character(id))
  lid <- levels(id)
  nmat <- nmat[lid, lid]

  return(list("map" = map, "nmat" = nmat, "data" = cbind(x, "id" = id)))
}


## Function to create a neighbormatrix from
## a list() of polygons or objects of class
## "SpatialPolygons".
neighbormatrix <- function(x, type = "dist", scale = NULL, k = 1) {
  require("maptools"); require("spdep")
  nx <- NULL
  ox <- x
  if(is.list(x)) {
    nx <- names(x)
    x <- bnd2sp(x)
  }
  if(is.null(type) || type == "delaunay" || type == "knear") {
    neighbors <- if(is.null(type)) {
      poly2nb(x)
    } else {
      if(type == "delaunay")
        tri2nb(coordinates(x))
      else
        knn2nb(knearneigh(coordinates(x), k = k, longlat = TRUE, RANN = FALSE), sym = FALSE)
    }
    adjmat <- nb2mat(neighbors, style = "B", zero.policy = TRUE)
    dims <- dim(adjmat)
    adjmat <- matrix(adjmat)
    dim(adjmat) <- dims
    if(!is.null(scale)) {
      require("fields")
      coords <- coordinates(x)
      codist <- fields::rdist(coords, coords)
      adjmat[(codist > (diff(range(codist)) * scale))] <- 0
    }
  }
  if(type == "dist") {
    require("fields")
    coords <- coordinates(x)
    codist <- fields::rdist(coords, coords)
    if(is.null(scale))
      scale <- 0.2
    adjmat <- 1 * (codist < (diff(range(codist)) * scale))
  }
  if(type == "boundary") {
    adjmat <- bnd2gra(ox)
    class(adjmat) <- "matrix"
    adjmat[adjmat != 0] <- 1
    diag(adjmat) <- 0
    nx <- NULL
  }
  if(!isSymmetric(adjmat)) {
    n <- nrow(adjmat)
    for(i in 1:n) {
      for(j in 1:n) {
        if(adjmat[i, j] > 0 & adjmat[j, i] < 1)
          adjmat[j, i] <- 1
        if(adjmat[i, j] < 1 & adjmat[j, i] > 0)
          adjmat[i, j] <- 1
      }
    }
  }
  if(is.null(nx))
    nx <- try(slot(x, "data")$NAME, silent = TRUE)
  if(is.null(nx))
    nx <- try(slot(x, "data")$ID, silent = TRUE)
  if(!is.null(nx) && class(nx) != "try-error") {
    rownames(adjmat) <- nx
    colnames(adjmat) <- nx
  }
  attr(adjmat, "coords") <- coordinates(x)
  adjmat
}


## Function to plot the neighborhood relationship
## given some polygon list() map x.
plotneighbors <- function(x, type = "dist", scale = NULL, k = 1,
  n.lwd = 1, n.col = "black", n.lty = 1, add = FALSE, ...) {
  adjmat <- neighbormatrix(x, type, scale, k)
  coords <- attr(adjmat, "coords")
  if(!add)
    plotmap(x, ...)
  args <- list(...)
  if(is.null(args$names)) points(coords)
  id <- 1:ncol(adjmat)
  for(i in 1:nrow(adjmat)) {
    neighbors <- id[adjmat[i, ] > 0]
    if(length(neighbors)) {
      xy.neighbors <- coords[c(i, neighbors), ]
      for(j in 2:nrow(xy.neighbors))
        lines(xy.neighbors[c(1, j), ], lwd = n.lwd, lty = n.lty, col = n.col)
    }
  }
  invisible(NULL)
}

