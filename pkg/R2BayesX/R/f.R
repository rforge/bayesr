f <-
function(x, z = NA, bs = "ps", by = NA, ...)
{
  rval <- list()
  rval$call <- call <- match.call()
  rval$term <- as.character(call$x)
  z <- deparse(substitute(z), backtick = TRUE, width.cutoff = 500)
  if(z != "NA") {
    rval$term <- c(rval$term, z)
    rval$dim <- 2L
  } else dim <- 1L
  rval$label <- paste("f(", paste(rval$term, sep = "", collapse = ","), ")", sep = "")
  rval$by <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
  xt <- list(...)
  rval$p.order <- rep(0, 2L)
  if(is.null(xt$degree))
    rval$p.order[1L] <- 3L
  else 
    rval$p.order[1L] <- xt$degree
  if(is.null(xt$order))
    rval$p.order[2L] <- 2L
  else
    rval$p.order[2L] <- xt$order
  if(is.null(xt$knots))
    rval$bs.dim <- -1L
  else 
    rval$bs.dim <- xt$knots + 2 * rval$p.order[1L] - 1
  xt[c("degree", "order", "knots", "nrknots")] <- NULL
  if(!length(xt))
    rval$xt <- NULL
  else
    rval$xt <- xt
  class(rval) <- paste(bs, ".smooth.spec", sep = "")

  return(rval) 
}

