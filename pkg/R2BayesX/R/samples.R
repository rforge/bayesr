samples <-
function(object, model = NULL, term = NULL, ...)
{
  if(is.null(term))
    term <- 1L
  object <- get.model(object, model)
  intcpt <- FALSE
  if(is.character(term) && !is.na(pmatch(term, "fixed-samples")))
    intcpt <- TRUE
  if(length(object) > 1L) {
    rval <- list()
    for(i in 1L:length(object)) {
      if(!intcpt) {
        rval[[i]] <- list(samples = attr(object[[i]]$effects[[term]], "sample"),
          variance.samples = attr(object[[i]]$effects[[term]], "variance.sample"))
      } else {
        rval[[i]] <- attr(object[[i]]$fixed.effects, "sample")
      }
    }
  } else {
    if(!intcpt) {
      rval <- list(samples = attr(object[[1L]]$effects[[term]], "sample"),
        variance.samples = attr(object[[1L]]$effects[[term]], "variance.sample"))
    } else {
      rval <- attr(object[[1L]]$fixed.effects, "sample")
    }
  }
  if(is.null(rval) || length(rval) < 1L)
    rval <- NA

  return(rval)
}

