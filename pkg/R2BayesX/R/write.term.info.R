write.term.info <-
function(file, terms, data, object = NULL, contrasts.arg = NULL, 
  xlev = NULL, intcpt = TRUE, rdafile)
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
          m <- model.matrix(as.formula(paste("~ -1 +", terms[k])), data, contrasts.arg, xlev)
          realname <- colnames(m)
          fn <- rmf(realname)
          realname <- paste("c(\'", paste(realname, sep = "", collapse = "\', \'"), "\')", sep = "")
          fnv <- "c("
          nf <- length(fn)
          if(nf > 1L)
            for(i in 1L:(nf - 1L))
              fnv <- paste(fnv, "\'", fn[i], "\',", sep = "")
          fnv <- paste(fnv, "\'", fn[nf], "\')", sep = "")
          xl <- paste(levels(x), collapse = "\',\'")
          xl <- paste("c(\'", xl , "\')", sep = "")
          info <- paste("list(term=\'", terms[k], "\',pos=" , k, 
            ",isFactor=TRUE", ",names=", fnv, ",levels=", xl, ",realname=", realname, ")", sep = "")
        } else { 
          info <- paste("list(term=\'", rmf(terms[k]), "\',pos=", k, ",isFactor=FALSE, realname=", 
            paste("\'", terms[k], "\'", sep = ""),")", sep = "")
        }
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
      f <- object$oformula
      save(f, file = rdafile)
      info <- paste("list(formula=\'", as.expression(f), "\',", sep = "")
      if(!is.null(object$hlevel))
        object$method <- "HMCMC"
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

