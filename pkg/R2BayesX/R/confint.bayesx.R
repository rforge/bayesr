confint.bayesx <- function(object, parm, level = 0.95, model = NULL, ...)
{
  object <- get.model(object, model)
  n <- length(object)
  if(n > 1L) {
    cf <- coef(object)
    medf <- edf(object)
    rval <- list()
    for(i in 1L:n) {
      if(!is.null(cf[[i]]))
        rval[[i]] <- confbayesx(cf[[i]], parm, level, medf[i], ...)
      else
        rval[[i]] <- NA
    }
    names(rval) <- names(object)
  } else {
    class(object) <- "bayesx"
    cf <- coef(object)
    medf <- edf(object)
    if(!is.null(cf))
      rval <- confbayesx(cf, parm, level, medf, ...)
    else
      rval <- NA
  }

  return(rval)
}


confbayesx <- function(x, parm, level = 0.95, edf, ...)
{
  cn <- colnames(x)
  pnames <- rownames(x)
  if(missing(parm)) 
    parm <- pnames
  else 
    if(is.numeric(parm)) 
      parm <- pnames[parm]
  id <- c(1L:length(pnames))[pnames %in% parm]
  ci <- NULL
  if(any(grepl("t value", cn))) {
    cf <- x[,1L]
    names(cf) <- pnames
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    fac <- qt(a, edf)
    pct <- paste(format(100 * a, trim = TRUE, 
      scientific = FALSE, digits = 3L), "%", sep = "")
    ci <- array(NA, dim = c(length(parm), 2L), 
      dimnames = list(parm, pct))
    ses <- x[id,2L]
    ci[] <- cf[parm] + ses %o% fac
  } else {
    samples <- attr(x, "sample")
    if(!is.null(samples)) {
      colnames(samples) <- rownames(x)
      samples <- samples[,id]
      if(!is.matrix(samples)) {
        samples <- matrix(samples, ncol = 1L)
        colnames(samples) <- parm
      }
      a <- (1 - level)/2
      a <- c(a, 1 - a)
      ci <- t(apply(samples, 2L, quantile, probs = a))
    }
  }

  return(ci)
}
