construct.shrw <- function(object, dir, prg, data, what)
{
  term <- object$term
  term <- paste(term, "(", what, sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA") {
    if(!missing(data) && !is.character(data)) {
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
