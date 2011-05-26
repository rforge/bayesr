write.term.info <-
function(file, terms, data, object = NULL, contrasts.arg = NULL, xlev = NULL)
{
  nt <- length(terms)
  if(nt > 0L) {
    for(k in 1L:nt) {
      if(is.sm(terms[k])) {
        te <- eval(parse(text = terms[k]))
        fby <- FALSE
        if(te$by != "NA") {
          if(!is.character(data)) {
            by <- data[[te$by]]
            if(is.factor(by)) {
              fby <- TRUE
              fnv <- paste("c(", paste("\'", te$by, levels(by), "\'", sep = "", collapse = ","), 
                ")", sep = "")
            }
          }
        }
        if(fby) {
          te$label <- gsub(")", paste(",by=", te$by, ")", sep = ""), te$label)
          info <- paste("list(term=\'", te$label, "\',pos=", k, ",by=\'", te$by,
            "\',isFactor=FALSE", ",isFactorBy=", fby, ",isFactorByNames=", fnv, ")", sep = "")
        } else {
          info <- paste("list(term=\'", te$label, "\',pos=", k, ",by=\'", te$by,
            "\',isFactor=FALSE", ",isFactorBy=", fby, ")", sep = "")
        }
      } else {
        sp <- FALSE
        if(grepl(":",terms[k]))
          sp <- TRUE
        if(!is.character(data) && !sp)
          x <- data[[terms[k]]]
        if(is.factor(x) && !sp) {
          m <- model.matrix(as.formula(paste("~", terms[k])), data, contrasts.arg, xlev)
          fn <- colnames(m)[2L:ncol(m)]
          fnv <- "c("
          nf <- length(fn)
          if(nf > 1L)
            for(i in 1L:(nf - 1L))
              fnv <- paste(fnv, "\'", fn[i], "\',", sep = "")
          fnv <- paste(fnv, "\'", fn[nf], "\')", sep = "")
          info <- paste("list(term=\'", terms[k], "\',pos=" , k, 
            ",isFactor=TRUE", ",names=",rmf(fnv), ")", sep = "")
        } else info <- paste("list(term=\'", terms[k], "\',pos=", k, ",isFactor=FALSE)", sep = "")
      }
      info <- paste(info,"\n")
      if(k < 2L)
        cat(info, file = file, append = FALSE)
      else
        cat(info, file = file, append = TRUE)
    }
    if(!is.null(object)) {
      if(!is.null(object$YLevels)) {
        YLevels <- paste(object$YLevels, collapse = "\\',\\'")
        YLevels <- paste("c(\\'", YLevels, "\\')", sep = "")
      }
      if(!is.null(object$nYLevels)) {
        nYLevels <- paste(object$nYLevels, collapse = "\\',\\'")
        nYLevels <- paste("c(\\'", nYLevels, "\\')", sep = "")
      }
      if(!is.null(object$order)) {
        ooo <- paste(object$order, collapse = ",")
        ooo <- paste("\'c(", ooo, ")\'", sep = "")
      }
      f <- as.character(object$oformula)
      f <- paste(f[2L], f[1L], f[3L])
      info <- paste("list(formula=\'", f, "\',", sep = "")
      info <- paste(info, "method=\'", object$method, "\',", sep = "")
      info <- paste(info, "family=\'", object$family, "\',", sep = "")
      info <- paste(info, "iterations=\'", object$iterations, "\',", sep = "")
      info <- paste(info, "step=\'", object$step, "\',", sep = "")
      if(!is.null(object$YLevels))
        info <- paste(info, "YLevels=\'", YLevels, "\',", sep = "")
      if(!is.null(object$nYLevels))
        info <- paste(info, "nYLevels=\'", nYLevels, "\',", sep = "")
      if(!is.null(object$order))
        info <- paste(info, "order=", ooo, ",", sep = "")
      info <- paste(info, "model.name=\'", object$model.name, "\')\n", sep = "")
      cat(info, file = file, append = TRUE)
    }
  }
  
  return(invisible(NULL))
}

