write.prg.setup <-
function(response, object, prg.file, data.file, thismodel, terms.specs)
{
  add.terms <- terms.specs$terms
  bt <- NULL
  vars <- terms.specs$vars
  if(!is.null(object$hlevel))
    if("(Intercept)" %in% vars)
      bt <- c(bt,"const")
  if(!is.null(add.terms)) {
    for(k in 1L:length(add.terms)) {
      st <- eval(parse(text = add.terms[k]))
      stclass <- class(st)
      if(stclass %in% c("tp.smooth.spec", "cs.smooth.spec")) {
        txt <- paste("constructor class \"", stclass, 
          "\" not supported by bayesx(), using bs = \"ps\" in term ", 
          st$label, "!", sep = "")
        warning(txt, call. = FALSE)
        class(st) <- "ps.smooth.spec"
      }
      bttmp <- bayesx.construct(st,object$outfile, prg.file, object$data)
      if(!is.null(object$hvariables))
        for(k in 1L:length(object$hvariables))
          if(bttmp == paste(object$hvariables[k], "(random)", sep = ""))
            bttmp <- paste(object$hvariables[k], "(hrandom)", sep = "")
      bt <- c(bt, bttmp)	
    }
  }
  bt <- c(bt, terms.specs$linvars)
  if("ModelOffset" %in% vars)
    bt <- c(bt, "ModelOffset(offset)")
  bt <- paste(bt, collapse = " + ")
  if("ModelWeights" %in% vars)
    bt <- paste(bt, "weight", "ModelWeights")
  control.values <- object[attr(object, "co.id")]
  predict <- control.values$predict
  control.values$predict <- NULL
  control.names <- names(control.values)
  if(is.null(object$hlevel))
    fullformula <- paste("b.regress ", response, " = ", bt, ",", sep = "")
  else
    fullformula <- paste("b.hregress ", response, " = ", bt, ",", sep = "")
  for(i in 1L:length(control.values))
    fullformula <- paste(fullformula, " ", control.names[i], "=", control.values[[i]], sep = "")
  dset <- paste("d", object$hlevel, sep = "")
  if(!is.null(object$hlevel)) {
    hlevel <- object$hlevel
    if(hlevel > 2L)
      hlevel <- 2L
    fullformula <- paste(fullformula, " hlevel=", hlevel, sep = "")
  }
  if(!is.null(predict))
    if(predict) {
      if(is.null(object$hlevel))
        fullformula <- paste(fullformula, "predict")
      else
        fullformula <- paste(fullformula, "predict=full")
    }
  fullformula <- paste(fullformula, "using", dset)
  cat("dataset", dset, "\n", file = prg.file, append = TRUE)
  cat(paste(dset, ".infile using ", data.file, "\n", sep = ""), file = prg.file, append = TRUE)
  cat(paste("b.outfile = ", object$outfile, "/", object$model.name, thismodel, "\n", sep = ""), 
    file = prg.file, append = TRUE)
  cat(fullformula, "\n", file = prg.file, append = TRUE)
  if(object$method == "MCMC" && object$first)
    cat("b.getsample\n", file = prg.file, append = TRUE)
  bayesx.prg <- paste(paste(readLines(paste(object$outfile, "/", prg.file, sep=""), n = -1L), 
    collapse = " \n"), " \n", sep="")

  ## only for hierarchical
  if(!is.null(object$hlevel)) {
    bayesx.prg <- gsub("psplinerw1", "pspline", bayesx.prg, fixed = TRUE)
    bayesx.prg <- gsub("psplinerw2", "pspline", bayesx.prg, fixed = TRUE)
    bayesx.prg <- gsub("maxint=150", "", bayesx.prg, fixed = TRUE)
    bayesx.prg <- gsub("aresp=1", "", bayesx.prg, fixed = TRUE)
    bayesx.prg <- gsub("bresp=0.005", "", bayesx.prg, fixed = TRUE)
    cat(bayesx.prg, file = prg.file)
  }

  return(bayesx.prg)
}

