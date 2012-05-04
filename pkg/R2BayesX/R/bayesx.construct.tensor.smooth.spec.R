bayesx.construct.tensor.smooth.spec <- bayesx.construct.t2.smooth.spec <-
function(object, dir, prg, data)
{
  by <- object$term[2L]
  term <- object$term[1L]
  object <- object$margin[[1L]]
  object$term <- term
  object$by <- by
  term <- bayesx.construct(object, dir, prg, data)
  term <- gsub("psplinerw2", "pspline2dimrw2", term)
  term <- gsub("psplinerw2", "pspline2dimrw1", term)

  return(term)
}

