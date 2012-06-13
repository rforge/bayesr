predict.bayesx <-
function(object, newdata, model = NULL, term = NULL, ...) 
{
  if(!missing(newdata))
    stop("out of sample prediction currently not available")

  return(fitted.bayesx(object, model, term, ...))
}

