sx <- function(x, z = NULL, bs = "ps", by = NULL, ...)
{
  call <- match.call()
  k <- -1
  m <- NA
  xt <- list(...)
  if(!is.null(xt$knots))
    xt$nrknots <- xt$knots
  if(bs %in% c("ps", "te", "psplinerw1", "psplinerw2", "pspline", "pspline2dimrw2")) {
    if(!is.null(xt$degree))
      m <- xt$degree
    if(!is.null(xt$order)) {
      if(is.na(m))
        m <- c(3L, xt$order)
      else
        m <- c(m[1L], xt$order)
    }
    if(length(m) < 2L && is.na(m))
      m <- c(3L, 2L)
    if(is.null(xt$order) && length(m) < 2L)
      m <- c(m, 2L)
    if(!is.null(xt$nrknots))
      k <- xt$nrknots + m[1L] - 1L
    else {
      if(bs %in% c("ps", "psplinerw1", "psplinerw2", "pspline"))
        k <- 20L + m[1L] + 1L
      else
        k <- 5L + m[1L] + 1L
    } 
  }
  if(bs %in% c("gs", "geospline")) {
    if(length(m) < 2L && is.na(m))
      m <- c(3L, 1L)
    else {
      if(length(m) < 2L)
        m <- c(m, 1L)
      else
        m[2L] <- 1L
    }
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
  if(is.null(by))
    rval$by <- "NA"
  else
    rval$by <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
  rval$label <- paste("sx(", paste(term, collapse = ",", sep = ""), ")", sep = "")

  return(rval)
}

