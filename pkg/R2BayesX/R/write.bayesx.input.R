write.bayesx.input <- function(object)
{
  if(is.null(object) || missing(object))
    stop("nothing to do!")
  if(class(object) != "bayesx.input")
    stop("object must be of class bayesx.input")
  data.file <- NULL
  if(is.null(object$data))
    object$data <- environment(object$formula)
  else {
    if(is.character(object$data))
      data.file <- object$data
  }
  if(is.null(object$outfile)) {
    if(.Platform$OS.type == "windows") {
      tmp <- splitme(tempdir())
      for(i in 1L:length(tmp))
        if(tmp[i] == "\\")
          tmp[i] <- "/"
      object$outfile <- paste(resplit(tmp), "/bayesx", sep = "")
    } else object$outfile <- paste(tempdir(), "/bayesx", sep = "")
  }
  if(!is.character(object$outfile))
    stop("argument outfile must be a character!")
  if(is.null(object$hlevel))
    is.h <- FALSE
  else {
    if(object$hlevel > 1L)
      is.h <- TRUE
    else
      is.h <- FALSE
  }
  if(!file.exists(object$outfile))
    dir.create(object$outfile, showWarnings = FALSE)
  else {
    if(!is.h && .Platform$OS.type == "unix") {
      wok <- TRUE
      count <- 0L
      while(wok) {
        count <- count + 1L
        if(!file.exists(paste(object$outfile, count, sep = ""))) {
          object$outfile <- paste(object$outfile, count, sep = "")
          wok <- FALSE
        }
      }
      dir.create(object$outfile, showWarnings = FALSE)
    }
  }
  wd <- as.character(getwd())
  setwd(object$outfile)
  prg.file <- paste(object$model.name, ".input.prg", sep = "")
  thismodel <- NULL
  if(file.exists(prg.file) && !is.h) {
    wok <- TRUE
    thismodel <- 0L
    while(wok) {
      thismodel <- thismodel + 1L
      prg.file <- paste(object$model.name, ".input", thismodel, ".prg", sep = "")
      if(!file.exists(prg.file))
        wok <- FALSE
    }
  }
  if(!is.h) {
    cat(paste("% usefile ", object$outfile, "/", prg.file, "\n", sep = ""), 
      file = prg.file, append = FALSE)
    cat(paste("logopen using ", object$outfile, "/", prg.file, ".log", "\n", sep = ""),
      file = prg.file, append = TRUE)
    if(object$method == "MCMC") {
      if(is.null(object$mcmcreg) && object$first)
        cat("bayesreg b\n", file = prg.file, append = TRUE)
      if(!is.null(object$mcmcreg) && object$first)
        cat("mcmcreg b\n", file = prg.file, append = TRUE)
    }
    if(object$method == "REML")
      cat("remlreg b\n", file = prg.file, append = TRUE)
    if(object$method == "STEP")
      cat("stepwisereg b\n", file = prg.file, append = TRUE)    
  }

  ## hierarchical models set here 
  if(!is.null(object$h.random))
    for(k in 1:length(object$h.random))
      write.bayesx.input(object$h.random[[k]])
  if(is.null(data.file)) {
    tf <- as.character(object$formula)
    if(grepl(tf[2L], tf[3L], fixed = TRUE))
      dat <- model.matrix(as.formula(paste(tf[1L], tf[3L])), object$data)
    else
      dat <- model.matrix(object$formula, object$data)
    nd <- rmf(names(object$data))
    names(object$data) <- nd
    colnames(dat) <- rmf(colnames(dat))
    dropvars <- NULL
    for(sf in nd)
      if(is.factor(ff <- object$data[[sf]])) {
        lf <- paste(sf, rmf(levels(ff)), sep = "")
        vars <- colnames(dat)
        if(!all(lf %in% vars)) {
          mf <- lf[!lf %in% vars]
          lf <- strsplit(mf, sf)
          dropvars <- c(dropvars, mf)
          for(k in 1:length(lf)) {
            dat <- cbind(dat, (lf[[k]][length(lf[[k]])] == ff) * 1)
            colnames(dat)[ncol(dat)] <- mf[k]
          }
        }
      }
    nc <- ncol(dat)
    if(!is.null(model.offset(object$data))) {
      nc <- ncol(dat)
      dat <- cbind(dat, model.offset(object$data))
      colnames(dat)[nc + 1L] <- "ModelOffset"
    }  
    if(!is.null(model.weights(object$data))) {
      nc <- ncol(dat)
      dat <- cbind(dat, model.weights(object$data))
      colnames(dat)[nc + 1L] <- "ModelWeights"
    }  
    dat <- cbind(as.vector(object$data[[object$Yn]]),dat)
    colnames(dat)[1L] <- object$Yn
    vars <- colnames(dat)
    nc <- ncol(dat)
    intcpt <- FALSE
    if("(Intercept)" %in% vars && nc > 1L) {
      dat <- matrix(dat[,!vars %in% "(Intercept)"], ncol = (nc - 1L))
      colnames(dat) <- vars[!vars %in% "(Intercept)"]
      intcpt <- TRUE
    }
    vars <- rmf(colnames(dat))
    for(i in 1L:length(vars)) {
      vars[i] <- gsub("(offset)", "ModelOffset", vars[i], fixed = TRUE)
      vars[i] <- gsub("(weights)", "ModelWeights", vars[i], fixed = TRUE)
    }
    colnames(dat) <- vars
    vars <- vars[vars!=object$Yn]
    if(intcpt)
      vars <- c("(Intercept)",vars)
    if(!is.null(object$hlevel)) {
      if(object$hlevel > 1L) {
        object$model.name <- paste(object$model.name, "_hlevel", 
          object$hlevel, "_RANDOM_EFFECTS", sep = "")
      } else {
        object$model.name <- paste(object$model.name, "_hlevel", 
          object$hlevel, "_MAIN_REGRESSION", sep = "")
      }
    }
    data.file <- paste(object$outfile, "/", object$model.name, ".data.raw", sep = "")
    if(!file.exists(data.file))
      write.table(dat, data.file, col.names = TRUE, row.names = FALSE, quote = FALSE)
    else {
      wok <- TRUE
      i <- 1L
      while(wok) {
        data.file <- paste(object$outfile, "/", object$model.name, ".data", i, ".raw", sep = "")
        i <- i + 1L
        if(!file.exists(data.file)) {
          write.table(dat, data.file, col.names = TRUE, row.names = FALSE, quote = FALSE)
          wok <- FALSE
        }
      }
    }
    terms <- attr(object$terms, "term.labels")
    infofile <- paste(object$outfile, "/", object$model.name, thismodel, sep = "")
    infofile <- paste(infofile, ".terms.info", sep = "")
    if(object$method != "STEP")
      write.term.info(infofile, terms, object$data, object)
  } else {
    vdf <- sub(" +$", "", readLines(data.file, n = 1L))
    vars <- strsplit(vdf, " ")[[1L]]
    ok <- FALSE
    if(length(vars) > 1L)
      if(all(attr(terms(object$formula), "term.labels") %in% vars))
        ok <- TRUE
    if(!ok)
      warning("variable names in specified data file do not match with formula variable names!")
  }

  st <- split.terms(attr(object$terms, "term.labels"), vars, object$data, dropvars)
  bayesx.prg <- write.prg.setup(object$Yn, object, prg.file, data.file,
    thismodel, st)
  if(!is.h)
    cat("logclose \n", file = prg.file, append = TRUE)
  setwd(wd)
  return(invisible(list(prg = bayesx.prg, prg.name = prg.file,
    model.name = object$model.name, file.dir = object$outfile)))
}


write.term.info <- function(file, terms, data, object = NULL, contrasts.arg = NULL, xlev = NULL)
{
  nt <- length(terms)
  if(nt > 0L) {
    for(k in 1L:nt) {
      if(is.sm(terms[k])) {
        te <- eval(parse(text = terms[k]))
        fby <- FALSE
        if(te$by != "NA") {
          by <- data[[te$by]]
          if(is.factor(by)) {
            fby <- TRUE
            fnv <- paste("c(", paste("\'", te$by, levels(by), "\'", sep = "", collapse = ","), 
              ")", sep = "")
          }
        }
        if(fby) {
          te$label <- gsub(")", paste(",by=", te$by, ")", sep = ""), te$label)
          info <- paste("list(term=\'", te$label, "\',pos=", k, ",by=\'", te$by,
            "\',isFactor=FALSE", ",isFactorBy=", fby, ",isFactorByNames=", fnv, ")", sep = "")
        } else {
          info <- paste("list(term=\'", te$label, "\',pos=", k, ",by=", te$by,
            ",isFactor=FALSE", ",isFactorBy=", fby, ")", sep = "")
        }
      } else {
        sp <- FALSE
        if(grepl(":",terms[k]))
          sp <- TRUE
        if(!sp)
          x <- data[[terms[k]]]
        if(is.factor(x) && !sp) {
          m <- model.matrix(as.formula(paste("~", terms[k])), data, contrasts.arg, xlev)
          fn <- colnames(m)[2L:ncol(m)]
          fnv <- "c("
          nf <- length(fn)
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
      f <- as.character(object$oformula)
      f <- paste(f[2L], f[1L], f[3L])
      info <- paste("list(formula=\'", f, "\',", sep = "")
      info <- paste(info, "method=\'", object$method, "\',", sep = "")
      info <- paste(info, "family=\'", object$family, "\',", sep = "")
      info <- paste(info, "iterations=\'", object$iterations, "\',", sep = "")
      info <- paste(info, "step=\'", object$step, "\',", sep = "")
      info <- paste(info, "model.name=\'", object$model.name, "\')\n", sep = "")
      cat(info, file = file, append = TRUE)
    }
  }
  
  return(invisible(NULL))
}


term.names <- function(x)
{
  n <- length(x)
  names <- vector("list",n)
  for(i in 1L:length(x)) {
    if(is.sm(x[i]))
      names[[i]] <- eval(parse(text = x[i]))$term
    else
      names[[i]] <- x[i]
  }
  
  return(names)
}


split.terms <- function(terms, vars, data, dropvars)
{
  if(length(terms) > 0L) {
    sm <- is.sm(terms)
    smt <- terms[sm]
    lv <- rmf(terms[!sm])
    if(length(smt) < 1L)
      smt <- NULL
    if(length(lv) < 1L)
      lv <- NULL
    linvars <- NULL
    for(i in 1L:length(lv)) {
      tmp <- gsub("(offset)", "ModelOffset", lv[i], fixed = TRUE)
      tmp <- gsub("(weights)", "ModelWeights", lv[i], fixed = TRUE)
      tmp2 <- data[[tmp]]
      if(is.factor(tmp2)) {
        lt <- rmf(levels(tmp2))
        tmp <- paste(tmp, lt, sep = "") 
        tmp <- tmp[tmp %in% vars]
        }
      linvars <- c(linvars, tmp)
    }
    if(!is.null(linvars) && !is.null(dropvars)) {
      for(drop in dropvars)
        linvars <- linvars[linvars != drop]
    }
  } else smt <- lv <- NULL
  if(!is.null(smt))
    smt <- na.omit(smt)
  if(!is.null(linvars))
    linvars <- na.omit(linvars)
  if(!is.null(vars))
    vars <- na.omit(vars)

  return(list(terms = smt, linvars = linvars, vars = vars))
}


is.sm <- function(x)
{
  n <- length(x)
  issm <- rep(FALSE, n)
  for(i in 1L:length(x)) {
    split <- splitme(x[i])
    if(length(split) > 2L) {
      if(resplit(split[1L:3L]) == "te(")
        issm[i] <- TRUE
      if(resplit(split[1L:2L]) %in% c("s(", "r("))
        issm[i] <- TRUE
    }
  }

  return(issm)
}


is.rt <- function(x)
{
  n <- length(x)
  isrt <- rep(FALSE, n)
  for(i in 1L:length(x)) {
    split <- splitme(x[i])
    if(length(split) > 2L) {
      if(resplit(split[1L:2L]) == c("r("))
        isrt[i] <- TRUE
    }
  }

  return(isrt)
}


## write to prg file
write.prg.setup <- function(response, object, prg.file, data.file, thismodel, terms.specs)
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


rmf <- function(x) 
{
  for(i in 1L:length(x))
    for(char in c("+", "-", "*", ":", "^", "/", " ")) 
      x[i] <- gsub(char, "_", x[i], fixed = TRUE)

  return(x)
}
