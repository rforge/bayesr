bayesx.construct.tensor.smooth.spec <- function(object, dir, prg, data)
{
  termo <- object$term
  object$p.order <- object$margin[[1L]]$p.order
  object$bs.dim <- (object$margin[[1L]]$bs.dim * 2L) + 1L
  if(is.na(object$p.order[1L]))
    object$p.order <- c(3L, 2L)
  if(object$bs.dim < 0L)
    object$bs.dim <- 22L
  else {
    if(length(object$p.order) > 1L) {
      if(object$p.order[2L] > 2L) {
        warning("order of the difference penalty not supported by BayesX, set to 2!")
        object$p.order <- c(object$p.order[1L], 2L)
      }
    }
  }
  nrknots <- object$bs.dim - object$p.order[1L] + 1L
  xt <- object$xt
  term <- paste(termo[1L], "*", termo[2L], "(pspline2dimrw", object$p.order[2L],
    ",nrknots=", nrknots, ",degree=", object$p.order[1L], sep = "")
  term <- paste(do.xt(term, object$margin[[1L]], NULL), ")", sep = "")
  if(object$by != "NA") {
    if(!is.character(data)) {
      by <- eval(parse(text = object$by), envir = data)
      if(is.factor(by))
        by <- paste(object$by, levels(by), sep = "")
      else
        by <- object$by
    } else by <- object$by
    term <- paste(by, "*", term, sep = "")
  }

  return(term)
}
