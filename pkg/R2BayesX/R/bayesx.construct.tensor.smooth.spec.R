bayesx.construct.tensor.smooth.spec <- bayesx.construct.t2.smooth.spec <-
function(object, dir, prg, data)
{
  termo <- object$term
  label <- object$label
  by <- object$by
  object <- object$margin[[1L]]
  if(class(object) != "ps.smooth.spec") {
    warning(paste("constructor class \"", class(object), 
      "\" not supported by bayesx(), using bs = \"ps\" in term ", label,"!", sep = ""), 
      call. = FALSE)
  }
  object$by <- by
  if(length(object$p.order) == 1L) 
    m <- rep(object$p.order, 2L)
  else 
    m <- object$p.order
  m[is.na(m)] <- 2L
  object$p.order <- m
  object$p.order[1L] <- object$p.order[1L] + 2L
  if(length(object$p.order) > 1L) {
    if(object$p.order[2L] > 2L) {
      warning("order of the difference penalty not supported by BayesX, set to 2!")
      object$p.order <- c(object$p.order[1L], 2L)
    }
  }
  if(object$bs.dim < 8L)
    object$bs.dim <- 8L
  nrknots <- object$bs.dim - object$p.order[1L] + 1L
  term <- paste(termo[1L], "*", termo[2L], "(pspline2dimrw", object$p.order[2L],
    ",nrknots=", nrknots, ",degree=", object$p.order[1L], sep = "")
  term <- paste(do.xt(term, object$xt, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}

