samples <-
function(object, model = NULL, term = NULL, acf = FALSE, ...)
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
          if(acf) {
            if(!is.null(tmp$coef.samples)) {              
              tmp$coef.samples.acf <- samplesacf(tmp$coef.samples, ...)
            }
            if(!is.null(tmp$variance.samples)) {
              tmp$variance.samples.acf <- samplesacf(tmp$variance.samples, ...)
            }
          }
          eval(parse(text = paste("rval[[i]]$'", tn[j], "' <- tmp", sep = "")))
        }
      } else {
        rval[[i]] <- list(coef.samples = attr(object[[i]]$effects[[term]], "sample"),
          variance.samples = attr(object[[i]]$effects[[term]], "variance.sample"))
        if(acf) {
          if(!is.null(rval[[i]]$coef.samples)) {
            rval[[i]]$coef.samples.acf <- samplesacf(rval[[i]]$coef.samples, ...)
          }
          if(!is.null(rval[[i]]$variance.samples)) {
            rval[[i]]$variance.samples.acf <- samplesacf(rval[[i]]$variance.samples, ...)
          }
        }
      }
    } else {
      if(acf) {
        if(!is.null(attr(object[[i]]$fixed.effects, "sample"))) {
          tmp <- list(coef.samples = attr(object[[i]]$fixed.effects, "sample"),
            coef.samples.acf = samplesacf(attr(object[[i]]$fixed.effects, "sample"), ...))
        } else tmp <- NULL
      } else  tmp <- attr(object[[i]]$fixed.effects, "sample")
      rval$'linear' <- tmp
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
  if(any(is.na(rval)))
    warning("samples are missing in object!")

  return(rval)
}


samplesacf <- function(x, ...) 
{
  if(is.matrix(x)) {
    rval <- NULL
    for(j in 1L:ncol(x))
      rval <- cbind(rval, stats::acf(stats::ts(x[,j]), plot = FALSE, ...)$acf)
    rval <- rval[2L:nrow(rval),]
    colnames(rval) <- colnames(x)
    rownames(rval) <- paste("lag-", 1L:nrow(rval), sep = "")
  } else {
    rval <- as.vector(stats::acf(stats::ts(x), plot = FALSE, ...)$acf)
    rval <- rval[2L:length(rval)]
    names(rval) <- paste("lag-", 1L:length(rval), sep = "")
  }

  return(rval)
}

