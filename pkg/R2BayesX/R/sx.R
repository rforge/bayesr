sx <- function(x, z = NULL, bs = "ps", by = NA, ...)
{
  by <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
  call <- match.call()
  k <- -1
  m <- NA
  xt <- list(...)
  if(by != "NA" && is.vector(by) && length(by) < 2L && (is.numeric(by) || is.integer(by))) {
    xt$b <- by
    by <- "NA"
  }
  if(!is.null(xt$knots))
    xt$nrknots <- xt$knots
  if(bs %in% c("ps", "te", "psplinerw1", "psplinerw2", "pspline",
    "pspline2dimrw2", "gs", "geospline")) {
    if(!is.null(xt$degree))
      m <- xt$degree
    if(!is.null(xt$order)) {
      if(is.na(m))
        m <- c(3L, xt$order)
      else
        m <- c(m[1L], xt$order)
    }
    if(length(m) < 2L && (bs %in% c("gs", "geospline")))
      m <- c(3L, 1L)
    if(length(m) < 2L && is.na(m))
      m <- c(3L, 2L)
    if(length(m) < 2L && is.na(m))
      m <- c(3L, 1L)    
    if(is.null(xt$order) && length(m) < 2L)
      m <- c(m, 2L)
    if(is.null(xt$order) && length(m) < 2L)
      m <- c(m, 1L)
    if(!is.null(xt$nrknots))
      k <- xt$nrknots + m[1L] - 1L
    else {
      if(bs %in% c("ps", "psplinerw1", "psplinerw2", "pspline"))
        k <- 20L + m[1L] - 1L
      else
        k <- 20L + m[1L] - 1L
    } 
    m[1L] <- m[1L] - 2L
  }
  if(bs %in% c("kr", "gk", "kriging", "geokriging")) {
    m <- c(1L, 1L)
    if(!is.null(xt$nrknots)) {
      k <- xt$nrknots
    } else k <- 20
  }
  if(!is.null(xt$map))
    xt$map.name <- as.character(call$map)
  xt[c("degree", "order", "knots", "nrknots")] <- NULL
  if(!length(xt))
    xt <- NULL
  term <- as.character(call$x)
  if(!is.null(call$z)) 
    term <- c(term, as.character(call$z))
  rval <- mgcv::s(x, z, k = k, bs = bs, m = m, xt = xt)
  rval$term <- term
  rval$by <- by
  rval$label <- paste("sx(", paste(term, collapse = ",", sep = ""), ")", sep = "")

  return(rval)
}

