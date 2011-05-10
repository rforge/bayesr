bayesx.construct.ps.smooth.spec <- function(object, dir, prg, data)
{
  if(is.na(object$p.order[1L]))
    object$p.order <- c(3L, 2L)
  if(length(object$p.order) < 2L)
    object$p.order <- c(object$p.order, 2L)
  if(object$bs.dim < 0L)
    object$bs.dim <- 10L
  else {
    if(length(object$p.order) > 1L) {
      if(object$p.order[2L] > 2L) {
        warning("order of the difference penalty not supported by BayesX, set to 2!")
        object$p.order <- c(object$p.order[1L], 2L)
      }
    }
  }
  nrknots <- object$bs.dim - object$p.order[1L] + 1L
  termo <- object$term
  xt <- object$xt
  term <- paste(termo, "(psplinerw", object$p.order[2L], ",nrknots=",
    nrknots, ",degree=", object$p.order[1L], sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}
