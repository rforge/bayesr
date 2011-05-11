AIC.bayesx <-
function(object, ..., k)
{
  args <- list(...)
  model <- args$model
  if(!is.null(args$print.names))
    print.names <- args$print.names
  else
    print.names <- FALSE 
  return(extract.model.diagnostic(object, model, "AIC", print.names))
}

