r <-
function(id, method = NULL, by = NA, xt = NULL, 
  data = NULL, weights = NULL, subset = NULL, 
  offset = NULL, na.action = na.fail, contrasts = NULL, 
  control = bayesx.control(...), ...)
{
  term <- deparse(substitute(id), backtick = TRUE, width.cutoff = 500L)
  by.var <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500L)
  ins <- formula <- NULL
  if(by.var == ".") 
    stop("by=. not allowed")
  if(term == ".") 
    stop("r(.) not yet supported.")
  call <- match.call()
  label <- paste("r(", term)
  if(!is.null(method)) {
    ins <- list()
    mlabel <- as.character(call[3L])
    split <- splitme(mlabel)
    if(split[1L] != "~")
      mlabel <- resplit(c("~", split))
    formula <- as.formula(paste(term,mlabel))
    label <- paste(label, ",", mlabel, collapse="")
    mf <- terms.formula(formula, specials=c("s", "te", "r"))
    mterms <- attr(mf, "term.labels")
    if(length(mterms) > 0L)
      for(k in 1L:length(mterms)) {
        if(is.sm(mterms[k]))
          ins[[k]] <- try(eval(parse(text = mterms[k])), silent = TRUE)
        else {
          ins[[k]] <-list(term = mterms[k], label = mterms[k])
          class(ins[[k]]) <- "lin.smooth.spec"
        }
      }
    }
  if(by.var != "NA")
    label <- paste(label, ",by=", by.var, collapse = "")
  label <- gsub(" ", "", paste(label, ")", sep = ""))
  rval <- list(term = term, label = label, by = by.var, xt = xt, 
    ins = ins, formula = formula, data = data, weights = weights, 
    subset = subset, offset = offset, na.action = na.action, 
    contrasts = contrasts, control = control)
  class(rval) <- "ra.smooth.spec"

  return(rval) 
}

