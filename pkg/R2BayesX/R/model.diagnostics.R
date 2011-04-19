BIC <- function(object, ...)
{
  UseMethod("BIC")
}


GCV <- function(object, ...)
{
  UseMethod("GCV")
}


DIC <- function(object, ...)
{
  UseMethod("DIC")
}


BIC.bayesx <- function(object, ...)
{
  args <- list(...)
  model <- args$model
  if(!is.null(args$print.names))
    print.names <- args$print.names
  else
    print.names <- FALSE 
  return(extract.model.diagnostic(object, model, "BIC", print.names))
}


GCV.bayesx <- function(object, ...)
{
  args <- list(...)
  model <- args$model
  if(!is.null(args$print.names))
    print.names <- args$print.names
  else
    print.names <- FALSE 
  return(extract.model.diagnostic(object, model, "GCV", print.names))
}


DIC.bayesx <- function(object, ...)
{
  args <- list(...)
  model <- args$model
  if(!is.null(args$print.names))
    print.names <- args$print.names
  else
    print.names <- FALSE 
  return(extract.model.diagnostic(object, model, "DIC", print.names))
}


AIC.bayesx <- function(object, ..., k)
{
  args <- list(...)
  model <- args$model
  if(!is.null(args$print.names))
    print.names <- args$print.names
  else
    print.names <- FALSE 
  return(extract.model.diagnostic(object, model, "AIC", print.names))
}


logLik.bayesx <- function(object, ...)
{
  args <- list(...)
  model <- args$model
  if(!is.null(args$print.names))
    print.names <- args$print.names
  else
    print.names <- FALSE 
  return(extract.model.diagnostic(object, model, "logLik", print.names))
}


edf <- function(object, model = NULL, print.names = FALSE)
{
  return(extract.model.diagnostic(object, model, "df", print.names))
}


pd <- function(object, model = NULL, print.names = FALSE)
{
  return(extract.model.diagnostic(object, model, "pd", print.names))
}


extract.model.diagnostic <- function(object, model, what, print.names = FALSE)
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
