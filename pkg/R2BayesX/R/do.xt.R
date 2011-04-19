do.xt <- function(term, object, not = NULL)
{
  if(!is.null(object$xt)) {
    names.xt <- names(object$xt)
    if(is.null(not))
      not <- "not"
    for(name in names.xt)
      if(!name %in% not) {
        if(name %in% c("full", "catspecific", "center", "derivative", "nofixed"))
          term <- paste(term, ",", name, sep = "")
        else
          term <- paste(term, ",", name, "=", object$xt[name], sep = "")
      }
  }

  return(term)
}
