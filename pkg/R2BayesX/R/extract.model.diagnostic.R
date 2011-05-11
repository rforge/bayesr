extract.model.diagnostic <-
function(object, model, what, print.names = FALSE)
{
  object <- get.model(object, model)
  dg <- NULL
  for(i in 1L:length(object)) {
    tmp <- eval(parse(text = paste("object[[i]]$model.fit$", what, sep = "")))
    if(is.null(tmp))
      tmp <- NA
    dg <- c(dg, tmp)
  }
  if(!is.null(dg) && print.names)
    names(dg) <- names(object)

  return(dg)
}

