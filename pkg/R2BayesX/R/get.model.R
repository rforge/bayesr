get.model <- function(x, model)
{
  if(!is.null(model)) {
    if(is.character(model)) {
      if(all(is.na(model <- pmatch(model, names(x)))))
        stop("argument model is specified wrong!")
    } else {
      if(max(model) > length(x) || is.na(model) || min(model) < 1) 
        stop("argument model is specified wrong!")
    }
    x <- x[model]
  }

  return(x)
}
