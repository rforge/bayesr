samples <-
function(object, model = NULL, term = NULL, ...)
{
  if(is.null(term))
    term <- 1L
  object <- get.model(object, model)
  intcpt <- FALSE
  if(is.character(term) && !is.na(pmatch(term, "linear-samples")))
    intcpt <- TRUE
  rval <- list()
  k <- length(object)
  mn <- rep("model", length.out = k)
  for(i in 1L:k) {
    if(!intcpt) {
      if(length(term) > 1L) {
        rval[[i]] <- list()
        tn <- names(object[[i]]$effects)
        if(is.null(tn))
          tn <- paste(term)
        for(j in 1:length(term)) {
          tmp <- list(coef.samples = attr(object[[i]]$effects[[term[j]]], "sample"),
            variance.samples = attr(object[[i]]$effects[[term[j]]], "variance.sample"))
          eval(parse(text = paste("rval[[i]]$'", tn[j], "' <- tmp", sep = "")))
        }
      } else {
        rval[[i]] <- list(coef.samples = attr(object[[i]]$effects[[term]], "sample"),
          variance.samples = attr(object[[i]]$effects[[term]], "variance.sample"))
      }
    } else {
      rval[[i]] <- attr(object[[i]]$fixed.effects, "sample")
    }
  if(!is.null(object[[i]]$bayesx.setup$model.name))
    mn[i] <- object[[i]]$bayesx.setup$model.name
  }
  mn[duplicated(mn)] <- paste(mn[duplicated(mn)], 1:length(mn[duplicated(mn)]) + 1L, sep = "")
  names(rval) <- mn
  if(length(rval) < 2L)
    rval <- rval[[1L]]
  if(is.null(rval) || length(rval) < 1L)
    rval <- NA

  return(rval)
}

