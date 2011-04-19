bayesx.construct.ra.smooth.spec <- function(object, dir, prg, data)
{
  term <- object$term
  if(is.null(object$ins))
    term <- paste(term, "(random", sep = "")
  else
    term <- paste(term, "(hrandom", sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA") {
    if(!missing(data) && !is.character(data)) {
      by <- eval(parse(text = object$by), envir = data)
      if(is.factor(by))
        by <- paste(object$by, levels(by), sep = "")
      else
        by <- object$by
    } else by <- object$by
  term <- paste(by, "*", term,sep = "")
  }

  return(term)
}
