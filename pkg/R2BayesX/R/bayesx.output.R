read.bayesx.output <- function(dir, model.name = NULL)
{
  if(is.null(dir))
    stop("no directory specified!")
  if(is.null(model.name))
    bm <- search.bayesx.models(dir)
  else
    bm <- model.name
  k <- length(bm)
  rval <- list()
  for(j in 1:k) {
    cmd <- paste("rval$", bm[j], "<-read.bayesx.model.output('", dir, "','", bm[j], "')", sep = "")
    eval(parse(text = cmd))
  }
  class(rval) <- "bayesx"

  return(rval)
}
	
	
search.bayesx.models <- function(dir)
{
  if(!file.exists(dir))
    stop("directory is not existing!")
  files <- list.files(dir)
  if(length(files) < 1L || is.null(files))
    stop(paste("no files in directory:",dir))
  ok <- FALSE
  model.names <- NULL
  if(length(i <- grep("bayesx.log", files, fixed = TRUE))) {
    logf <- readLines(paste(dir,"/",files[i],sep=""))
    outfiles <- grep("b.outfile = ", logf, value = TRUE)
    if(length(outfiles) > 0L)
      for(k in 1L:length(outfiles)) {
        of <- splitme(outfiles[k])
        go <- TRUE
        tk <- NULL
        for(j in length(of):1) {
          if(of[j] == "/" || of[j] == "\\")
            go <- FALSE
          if(go)
            tk <- c(tk,of[j])  
        }
        model.names <- c(model.names,resplit(tk[length(tk):1]))
      }
    collect <- NULL
    for(k in 1L:length(model.names))
      if(any(grepl(model.names[k], files)))
        collect <- c(collect, model.names[k])
    if(!is.null(collect)) {
      ok <- TRUE
      model.names <- collect
    }
  }

  ## specifiy possible search endings here
  search.endings <- c("_model_summary.tex","_graphics.prg","_stata.do",
    "_r.R","_FixedEffects.res","_FixedEffects1.res","_FixedEffects2.res",
    "_predict.raw","_modelfit.raw","_LinearEffects.res","_LinearEffects1.res",
    "_LinearEffects2.res",".terms.info","_predict.res")

  model.names2 <- list()
  n <- 0L
  for(i in 1L:length(search.endings)) {
    if(length(ii <- grep(search.endings[i],files))) {
      n <- n + 1L
      fs <- files[ii]
      k <- length(fs)
      model.names2[[n]] <- rep("NA",k)
      for(j in 1L:length(fs))
        model.names2[[n]][j] <- strsplit(fs[j], search.endings[i], "")[[1L]]
      ok <- TRUE
    }
  }
  model.names <- c(model.names, unlist(model.names2))
  if(!ok)
    stop(paste("no BayesX output files found in:",dir))
  else {
    model.names <- unique(unlist(model.names))
    for(k in 1L:length(model.names)) {
      model.names[k] <- strsplit(model.names[k], "_MAIN_REGRESSION", "")[[1L]][1L]
      model.names[k] <- strsplit(model.names[k], "_RANDOM_EFFECTS", "")[[1L]][1L]
      if(length(grep("hlevel",model.names[k]))) {
        split <- strsplit(model.names[k],"_hlevel","")[[1L]]
        if(length(split) < 2L)
          split <- strsplit(model.names[k],".hlevel","")[[1L]]
        mn <- split[1L]
        hlevel <- split[2L]
        if(hlevel < 2L)
          model.names[k] <- paste(mn,"_hlevel",hlevel,"_MAIN_REGRESSION",sep="")
        else
          model.names[k] <- paste(mn,"_hlevel",hlevel,"_RANDOM_EFFECTS",sep="")
      }
    }
  }
		
  return(unique(unlist(model.names)))
}
	
	
splitme <- function(x)
{
  return(strsplit(x, "")[[1L]])
}
	
	
resplit <- function(x)
{
  if(!is.null(x))
    x <- paste(x, sep = "", collapse = "")
	
  return(x)
}

	
search.bayesx.tex <- function(x)
{
  rval <- list()
  fam <- any(grepl("Family:", x))
  if(fam) {
    fam <- grep("Family:", x, value = TRUE)
    fam <- splitme(fam)
    nfam<- length(fam)
    j <- nfam + 1L
    ifam <- NULL
    for(i in 1L:nfam) {
      if(fam[i] == ">")
        j <- i
      if(fam[i] == "\\")
        j <- nfam + 1L
      if(i > j)
        ifam <- c(ifam, fam[i])
    }
    rval$family <- resplit(ifam)
  }
  m <- any(grepl("BAYESREG", x))
  if(m)
    m <- "MCMC"
  else {
    m <- any(grepl("remlreg", x))
    if(m)
      m <- "REML"
    else {
      m <- any(grepl("STEPWISEREG", x))
      if(m)
        m <- "STEP"
      else
        m <- NULL
    }
  }
  rval$method <- m
  obs <- any(grepl("Number of observations", x))
  if(obs) {
    obs <- grep("Number of observations:", x, value = TRUE)
    obs <- splitme(obs)
    nobs <- length(obs)
    j <- nobs + 1L
    N <- NULL
    for(i in 1L:nobs) {
      if(obs[i] == "=" || obs[i] == ">")
        j <- i
      if(obs[i]=="\\")
        j <- nobs + 1L
      if(i > j)
        N <- c(N, obs[i])
    }
    rval$N <- as.integer(resplit(N))
    if(is.na(rval$N))
      rval$N <- NULL
  }
  iter <- any(grepl("Number of Iterations:", x))
  if(iter) {
    iter <- splitme(grep("Number of Iterations:", x, value = TRUE))
    ni <- length(iter)
		j <- ni + 1L
    igrep <- NULL
    for(i in 1L:ni) {
      if(iter[i] == ">")
        j <- i
      if(iter[i] == "\\")
        j <- ni + 1L
      if(i > j)
        igrep <- c(igrep, iter[i])
    }
    rval$iterations <- as.integer(resplit(igrep))
  }
  burn <- any(grepl("Burn in:", x))
  if(burn) {
    burn <- splitme(grep("Burn in:", x, value = TRUE))
    ni <- length(burn)
    j <- ni + 1L
    burngrep <- NULL
    for(i in 1L:ni) {
      if(burn[i] == ">")
        j <- i
      if(burn[i] == "\\")
        j <- ni + 1L
      if(i > j)
        burngrep <- c(burngrep,burn[i])
    }
    rval$burnin <- as.integer(resplit(burngrep))
  }
  thin <- any(grepl("Thinning Parameter:", x))
  if(thin) {
    thin <- splitme(grep("Thinning Parameter:", x, value = TRUE))
    ni <- length(thin)
    j <- ni + 1L
    thingrep <- NULL
   for(i in 1L:ni) {
     if(thin[i] == ">")
       j <- i
     if(i > j)
       thingrep <- c(thingrep, thin[i])
    }
    rval$thinning <- as.integer(resplit(thingrep))
  }
  dic <- any(grepl("DIC \\\\>", x))
  if(dic) {
    dic <- splitme(grep("DIC \\\\>", x, value = TRUE))
    ni <- length(dic)
    j <- ni + 1L
    dicgrep <- NULL
    for(i in 1L:ni) {
      if(dic[i] == ">")
        j <- i
      if(dic[i] == "\\")
        j <- ni + 1L
      if(i > j)
        dicgrep <- c(dicgrep, dic[i])
    }
    rval$DIC <- as.numeric(resplit(dicgrep))
  }
  pd <- any(grepl("pD \\\\>", x))
  if(pd) {
    pd <- splitme(grep("pD \\\\>", x, value = TRUE))
    ni <- length(pd)
    j <- ni + 1L
    pdgrep <- NULL
    for(i in 1L:ni) {
      if(pd[i] == ">")
        j <- i
      if(pd[i]=="\\")
        j <- ni + 1L
      if(i > j)
        pdgrep <- c(pdgrep, pd[i])
    }
    rval$pd <- as.numeric(resplit(pdgrep))
  }
  df <- any(grep("Degrees of freedom:", x))
  if(df) {
    df <- grep("Degrees of freedom:", x, value = TRUE)
    if(length(df) > 1L)
      df <- df[2L]
    df <- splitme(df)
    ok <- FALSE
    DF <- NULL
    for(i in 1L:length(df)) {
      if(ok && df[i] != "\\")
        DF <- c(DF,df[i])
      if(df[i] == ">")
        ok <- TRUE
    }
    rval$df <- as.numeric(resplit(DF))
  }
  aic <- any(grep("AIC:", x))
  if(aic) {
    aic <- grep("AIC:", x, value = TRUE)
    aic <- splitme(aic)
    ok <- FALSE
    AIC <- NULL
    for(i in 1L:length(aic)) {
      if(ok && aic[i] != "\\")
        AIC <- c(AIC, aic[i])
      if(aic[i] == ">")
        ok <- TRUE
    }
    rval$AIC <- as.numeric(resplit(AIC))
  }
  bic <- any(grep("BIC:", x))
  if(bic) {
    bic <- grep("BIC:", x, value = TRUE)
    bic <- splitme(bic)
    ok <- FALSE
    BIC <- NULL
    for(i in 1L:length(bic)) {
      if(ok && bic[i] != "\\")
        BIC <- c(BIC, bic[i])
      if(bic[i] == ">")
        ok <- TRUE
    }
    rval$BIC <- as.numeric(resplit(BIC))
  }
  gcv <- any(grep("GCV:", x))
  if(gcv) {
    gcv <- grep("GCV:", x, value = TRUE)
    gcv <- splitme(gcv)
    ok <- FALSE
    GCV <- NULL
    for(i in 1L:length(gcv)) {
      if(ok && gcv[i] != "\\")
        GCV <- c(GCV, gcv[i])
      if(gcv[i] == ">")
        ok <- TRUE
    }
    rval$GCV <- as.numeric(resplit(GCV))
  }
  final <- any(grepl("Final Predictor:", x))
  if(final) {
    final <- i <- grep("Final Predictor:", x)
    stepfiles <- NULL
    run <- TRUE
    while(run) {
      i <- i + 1L
      if(x[i] == "\\newpage ")
        run <- FALSE
      else
        stepfiles <- c(stepfiles, x[i])
      if(i == length(x))
        run <- FALSE
    }
    final <- splitme(grep("eta",stepfiles,value=TRUE))
    grepfinal <- NULL
    for(i in 1L:length(final)) {
      check <- final[i] != "$" && final[i] != "&" && final[i] != "\\"
      if(check) {
        ok <- final[i]
        if(ok == " ")
          ok <- "!"
        grepfinal <- c(grepfinal, ok)
      }
    }
    grepfinal <- resplit(grepfinal)
    grepfinal <- sub("cdot", "*", grepfinal)
    grepfinal <- sub("eta..", "eta!", grepfinal)
    grepfinal <- sub("..gamma_0", "!gamma0", grepfinal)
    grepfinal <- splitme(grepfinal)
    stepgrep <- NULL
    ok <- TRUE
    for(i in 1L:length(grepfinal)) {
      if(grepfinal[i] == "{")
        ok <- FALSE
      if(ok) {
        take <- grepfinal[i]
        if(take == "!")
          take <- " "
        if(take == "_")
          take <- NULL
        stepgrep <- c(stepgrep, take)
      }
      if(grepfinal[i] == "}")
        ok <- TRUE	
    }
    crit <- splitme(x[grep("Final Predictor:", x) + 8L])
    fc <- fcn <- NULL
    okn <- TRUE
    okc <- FALSE
    for(i in 1L:length(crit)) {
      if(okc && crit[i] != "\\")
        fc <- c(fc, crit[i])
      if(crit[i] == "=") {
        okn <- FALSE
        okc <- TRUE
      }
      if(okn && crit[i] != "\\")
        fcn <- c(fcn, crit[i])
    }
    fc <- as.numeric(resplit(fc))
    fcn <- resplit(fcn)
    fcn <- sub(" ", "", fcn)
    eval(parse(text=paste("rval$", fcn, "<-", fc, sep = "")))
    rval$step.final.model <- resplit(stepgrep)	
  }

  return(rval)
}	


read.bayesx.model.output <- function(dir, model.name)
{
  if(is.null(dir))
    stop("no directory specified!")
  if(is.null(model.name))
    stop("no model name specified!")
  files <- dir.files <- list.files(dir)
  if(!any(grep(model.name,files)))
    stop(paste("no model results existing for ", model.name, "!", sep = ""))
  else {
    rval <- list()
    files <- grep(model.name, files, value = TRUE, fixed = TRUE)
    filep <- grep(paste(model.name, ".", sep = ""), files, value = TRUE, fixed = TRUE)  
    files <- c(grep(paste(model.name, "_", sep = ""), files, value = TRUE, fixed = TRUE), filep)  

    ## search for data
    data <- NULL
    if(length(i <- grep(paste(model.name, ".data.raw", sep = ""), files)))
      data <- as.matrix(read.table(paste(dir, "/", files[i], sep = ""), header = TRUE))
    for(char in c("_predict.raw", "_predictmean.raw", "_predict.res"))
      if(length(i <- grep(char, files)))
        data <- cbind(data, as.matrix(read.table(paste(dir, "/", files[i], sep = ""), header = TRUE)))

    ## set response and predictor
    N <- NA
    response <- eta <- residuals <- NULL
    if(!is.null(data)) {
      data[data == "."] <- "NA"
      mode(data) <- "numeric"
      data <- as.data.frame(data)
      data$intnr <- NULL
      dn <- unique(names(data))
      data <- data[dn]
      if(is.null(response))
        response <- data[[1L]]
      if(is.factor(response)) 
        response <- as.numeric(as.character(response))
      N <- length(response)
      rval$fitted.values <- eta <- get.eta(data)
      rval$residuals <- response - eta
      rval$response <- response
    }

    ## get smooth and random effects
    rval <- c(rval, find.smooth.random(dir, files, data, response, eta, model.name))

    ## get fixed effects
    rval <- find.fixed.effects(dir, files, data, response, eta, model.name, rval)

    ## get scale estimate
    rval$variance <- get.scale(files, dir)

    ## search for other results
    model.results <- mf <- mf2 <- NULL
    method <- ""
    if(any(grep("deviance.raw", files))) {
      mf <- grep("deviance.raw", files, value = TRUE)
      mf <- read.table(paste(dir, "/", mf, sep = ""), header = TRUE)
      pd <- mf$unstandardized_deviance[length(mf$unstandardized_deviance) - 1L]
      DIC <- mf$unstandardized_deviance[length(mf$unstandardized_deviance)]
      mf <- list(DIC = DIC, pd = pd)
    }
    if(any(grep(".tex",files))) {
      sm <- readLines(paste(dir, "/", grep(".tex", files, value = TRUE), sep = ""))
      model.results <- search.bayesx.tex(sm)
      method <- model.results$method
    }
    if(any(grep("modelfit.raw",files))) {
      mf <- grep("modelfit.raw", files, value = TRUE)
      mf <- chacol(read.table(paste(dir, "/", mf, sep = ""), header = TRUE))
      mf <- mf2 <- as.list(mf)
    }
    if(!is.null(model.results)) {
      if(!is.null(mf2)) {
        n1 <- names(mf2)
        n2 <- names(model.results)
        for(i in 1L:length(mf2))
          if(!n1[i] %in% n2)
            eval(parse(text = paste("model.results$", n1[i], "<- mf2[[i]]", sep = "")))
      }
      mf <- model.results
    } else mf <- c(list(method = method, N = N), mf)
    if(mf$method == " ")
      mf$method <- "NA"
    if(!is.null(mf$logLik))
      mf$logLik <- mf$logLik/(-2)
    rval$model.fit <- mf

    ## reordering and naming
    if((info <- paste(model.name, ".terms.info", sep = "")) %in% dir.files)
      rval$effects <- term.reorder(rval$effects, paste(dir, "/", info, sep = ""))

    ## search for additional info
    rval$model.fit <- smi(info, rval$model.fit)

    ## reformate output
    rval <- bayesx.reformate(rval)

    return(rval)
  }
}


smi <- function(info.file, model.fit)
{
  info <- NULL
  if(file.exists(info.file))
    if(length(info <- readLines(info.file))) {
      if(is.null(model.fit))
        model.fit <- eval(parse(text = info[length(info)]))
      else {
        model.fit2 <- eval(parse(text = info[length(info)]))
        nmf2 <- names(model.fit2)
        nmf <- names(model.fit)
        model.fit[nmf[nmf %in% nmf2]] <- NULL
        model.fit <- c(model.fit, model.fit2)
      }
    }
  if(!is.null(model.fit$family)) {
    if(model.fit$family == " Gaussian" || model.fit$family == "Gaussian")
      model.fit$family <- "gaussian"
  }
  if(model.fit$method == "REML") {
    model.fit$step <- NULL
    model.fit$iterations <- NULL
  }

  return(model.fit)
}


df2m <- function(x)
{
  x$intnr <- NULL
  x <- as.matrix(x)
  rownames(x) <- 1L:nrow(x)

  return(x)
}


df2m2 <- function(x)
{
  x$intnr <- NULL
  x$paramnr <- NULL
  vars <- as.character(x$varname)
  x$varname <- NULL
  x <- as.matrix(x)
  vars[vars == "const"] <- "(Intercept)"
  rownames(x) <- vars

  return(x)
}


s4dim <- function(x)
{
  dim <- 1L
  if(ncol(x) == 10L)
    dim <- 2L
  if(ncol(x) == 11L)
    dim <- 2L
  if(ncol(x) == 16L)
    dim <- 2L

  return(dim)
}


s4class <- function(x)
{
  if(grepl("_random", x))
    cx <- "random.bayesx"
  if(grepl("_pspline", x))
    cx <- "sm.bayesx"
  if(grepl("_rw", x))
    cx <- "sm.bayesx"
  if(grepl("_spatial", x))
    cx <- "mrf.bayesx"
  if(grepl("_geospline", x))
    cx <- "geo.bayesx"
  if(grepl("_geokriging", x))
    cx <- "geo.bayesx"
  if(grepl("_logbaseline", x))
    cx <- "sm.bayesx"
  if(grepl("_kriging", x))
    cx <- "sm.bayesx"

  return(cx)
}


make.label <- function(cx, xnam, dimx, vx)
{
  if(cx == "random.bayesx")
    label <- paste("r(", xnam, sep = "")
  if(cx == "sm.bayesx") {
    if(dimx > 1L)
      label <- paste("s(", xnam[1L], ",", xnam[2L], sep = "")
    else
      label <- paste("s(", xnam, sep = "")
  }
  if(cx == "mrf.bayesx")
    label <- paste("s(", xnam, sep = "")
  if(cx == "geo.bayesx")
    label <- paste("s(", xnam[1], sep = "")
  if(is.null(vx))
    label <- paste(label, ")", sep = "")
  else
    label <- paste(label, "):", vx, sep = "")

  return(label)
}


blow.up.resid <- function(data, x, xnam, response, eta, dimx, cx)
{
  if(!is.null(data)) {
    x <- x[order(x[,1L]),]
    if(!is.matrix(x))
      x <- matrix(x, nrow = 1)
    id <- NULL
    if(cx == "geo.bayesx") {
      co <- x[,2L:3L]
      id <- x[,1L]
    } else co <- matrix(x[,1L:dimx], ncol = dimx)
    xtmp <- data[[xnam[1L]]]
    ind <- unique.id(xtmp)[order(xtmp)]
    x <- as.data.frame(x[ind,])
    x$pcat80 <- NULL
    x$pcat95 <- NULL
    x$pcat80_sim <- NULL
    x$pcat95_sim <- NULL
    x$pcat80_tot <- NULL
    x$pcat95_tot <- NULL
    x$pcat80tot_sim <- NULL
    x$pcat95tot_sim <- NULL
    for(k in 1L:length(xnam))
      eval(parse(text = paste("x$", xnam[k] , "<- NULL", sep = "")))
    x <- as.matrix(x)
    pres <- response - eta[,1L] + x[,1L]
    pres <- pres - mean(pres, na.rm = TRUE) ## mean(x[,1L], na.rm = TRUE)
    x <- cbind(co[ind,], pres, id[ind])
    if(ncol(x) < 3L)
      colnames(x) <- c("x.co", "partial.resids")
    else {
      if(!is.null(id))
        colnames(x) <- c("x.co", "y.co", "partial.resids", "id")
      else
        colnames(x) <- c("x.co", "y.co", "partial.resids")
    }
    rownames(x) <- 1L:nrow(x)
  }

  return(x)
}


get.eta <- function(data)
{
  etaspec <- c("eta","linpred","pmean_pred","pqu2p5_pred","pqu10_pred",
    "pmed_pred","pqu90_pred","pqu97p5_pred","pmean_mu","pqu2p5_mu",
    "pqu10_mu","pmed_mu","pqu90_mu","pqu97p5_mu","mu")
  nd <- names(data)
  eta <- enam <- NULL
  for(e in etaspec) 
    for(n in nd)
      if(!is.na(which <- pmatch(e, n))) {
        eta <- cbind(eta, as.matrix(data[n], mode = "numeric")) 
        enam <- c(enam, n)
      }
  if(!is.null(eta)) {
    colnames(eta) <- enam
    rownames(eta) <- 1L:nrow(eta)
    storage.mode(eta) <- "numeric"
  }

  return(eta)
}


find.smooth.random <- function(dir, files, data, response, eta, model.name)
{
  effects <- list()
  SmoothHyp <- RandomHyp <- NULL
  if(any((i <- grep(".res", files)))) {
    resfiles <- files[i]
    endings <- c("_predict.res", "FixedEffects", "LinearEffects", "scale.res", "_variance_", 
      "_var.res", ".raw", "_param.res", "_interact.res", "_df.res", "_lambda.res", 
      "_knots.raw", "_contour.res")
    for(res in endings)
      resfiles <- resfiles[!grepl(res, resfiles)]
    resfiles <- resfiles[!grepl("_lasso", resfiles)]
    resfiles <- resfiles[!grepl("_ridge", resfiles)]
    resfiles <- resfiles[!grepl("_nigmix", resfiles)]
    if(length(resfiles) > 0L)
      for(res in resfiles) {
        x <- df2m(read.table(paste(dir, "/", res, sep = ""), header = TRUE))
        dimx <- s4dim(x)
        if(sum(x[,(dimx + 1L):ncol(x)], na.rm = TRUE) != 0) {
          x <- x[order(x[,1L]),]
          cx <- s4class(res)
          dimx2 <- dimx
          if(cx == "geo.bayesx") {
            dimx2 <- dimx + 1L
            xnam <- colnames(x)[1L:dimx2]
            xnam2 <- xnam[1L]
          } else xnam <- xnam2 <- colnames(x)[1L:dimx2]
          xnam <- unlist(strsplit(xnam, "_c"))
          xnam2 <- xnam
          vx <- NULL
          if(grepl("_f_", res)) {
            res2 <- gsub(paste(model.name, "_", sep = ""), "", strsplit(res, "_f_")[[1L]])[1L]
            if(res2 != model.name)
              vx <- res2
          }
          colnames(x)[1L:dimx2] <- xnam
          rownames(x) <- 1L:nrow(x)
          labelx <- make.label(cx, xnam, dimx, vx)
          if(grepl("_spatialtotal.res", res))
            labelx <- paste(labelx, ":total", sep = "")
          attr(x, "partial.resids") <- blow.up.resid(data, x, xnam, response, eta, dimx, cx)
          attr(x, "specs") <- list(dim = dimx, term = xnam, 
            by = vx, label = labelx, is.factor = FALSE)
          ## search and set additional attributes
          nx <- length(xnam2)
          if(nx > 1L) {
            if(nx > 2L) {
              af1 <- grep(paste("_", xnam[1L], ".res", sep = ""), files, value = TRUE)
              af2 <- grep(paste("_", xnam2[1L], "_", sep = ""), files, value = TRUE)
            } else {
              af1 <- grep(paste("_", xnam[1L], "_", xnam[2L],".res", sep = ""), files, value = TRUE)
              af2 <- grep(paste("_", xnam2[1L], "_", xnam[2L], sep = ""), files, value = TRUE)
            }
          } else {
            af1 <- grep(paste("_", xnam[1L], ".res", sep = ""), files, value = TRUE)
            af2 <- grep(paste("_", xnam2[1L], "_", sep = ""), files, value = TRUE)
          }
          af <- c(af1, af2)
          if(!is.null(vx))
            af <- af[grepl(paste(vx, "_", sep = ""), af)]
          if(any(grep("_random", res, fixed = TRUE)))
            af <- grep("_random", af, fixed = TRUE, value = TRUE)
          if(any(grep("_spatial", res, fixed = TRUE)))
            af <- grep("_spatial", af, fixed = TRUE, value = TRUE)
          if(any(grep("_geokriging", res, fixed = TRUE)))
            af <- grep("_geokriging", af, fixed = TRUE, value = TRUE)
          if(length(af) > 0L) {
            if(length(varf <- grep("_var", af, value = TRUE))) {
              if(length(vf <- grep("_var.res", varf, value = TRUE))) {
                attr(x, "variance") <- df2m(read.table(paste(dir, "/", vf, sep = ""), 
                  header = TRUE))
                rownames(attr(x, "variance"))[1] <- labelx
                if(cx == "random.bayesx")
                  RandomHyp <- rbind(RandomHyp, attr(x, "variance"))
                else
                  SmoothHyp <- rbind(SmoothHyp, attr(x, "variance"))
              }
              if(length(vf <- grep("_variance_", varf, value = TRUE))) {
                if(length(vf2 <- vf[!grepl("sample", vf)])) {
                  attr(x, "variance") <- df2m(read.table(paste(dir, "/", vf2, sep = ""), 
                    header = TRUE))
                  rownames(attr(x, "variance"))[1] <- labelx
                  if(cx == "random.bayesx")
                    RandomHyp <- rbind(RandomHyp, attr(x, "variance"))
                  else
                    SmoothHyp <- rbind(SmoothHyp, attr(x, "variance"))
                }
                if(length(vf2 <- vf[grepl("sample", vf)]))
                  for(tf in vf2) {
                    attr(x, "variance.sample") <- df2m(read.table(paste(dir, "/", tf, sep = ""), 
                      header = TRUE))
                  }
              }
            }
            if(length(sf <- grep("_sample", af, value = TRUE)))
              if(length(sf <- sf[!grepl("_variance_", sf)]))
                for(tf in sf) {
                  attr(x, "sample") <- df2m(read.table(paste(dir, "/", tf, sep = ""), 
                    header = TRUE))
                }
            if(length(pf <- grep("_param", af, value = TRUE)))
              for(tf in pf)
                attr(x, "param") <- df2m(read.table(paste(dir, "/", tf, sep = ""), header = TRUE))
            if(length(kf <- grep("_knots", af, value = TRUE)))
              attr(x, "knots") <-  df2m(read.table(paste(dir, "/", kf, sep = ""), header = TRUE))
            if(length(cf <- grep("_contour", af, value = TRUE))) {
              attr(x, "contourprob") <-  df2m(read.table(paste(dir, "/", cf, sep = ""), 
                header = TRUE))
            }
          }
          class(x) <- cx
          eval(parse(text = paste("effects$\'", labelx, "\' <- x", sep = "")))
        }
      }      
  }

  return(list(effects = effects, smooth.hyp = mum(SmoothHyp), random.hyp = mum(RandomHyp)))
}


mum <- function(x)
{
  rn <- rownames(x)
  a <- duplicated(x[,1L], fromLast = FALSE)
  b <- duplicated(x[,1L], fromLast = TRUE)
  a[a != b] <- TRUE
  if(any(tot <- grepl(":total", rn))) {
    x <- x[!tot,]
    rn <- rn[!tot]
  }
  if(any(a)) {
    drn <- rn[a]
    nc <- nchar(drn)
    rmrn <- drn[nc < max(nc, na.rm = TRUE)]
    x <- x[!(rn %in% rmrn),]
  }
    
  return(x)
}


find.fixed.effects <- function(dir, files, data, response, eta, model.name, rval)
{
  FixedEffects <- NULL
  fefiles <- c(grep("LinearEffects", files, value = TRUE), 
    grep("FixedEffects", files, value = TRUE))
  lasso <- grep("_lasso_Effects", files, value = TRUE)
  ridge <- grep("_ridge_Effects", files, value = TRUE)
  nigmix <- grep("_nigmix_Effects", files, value = TRUE)
  fefiles <- c(fefiles, lasso, ridge, nigmix)
  if(length(fefiles)) {
    fsample <- fdf <- NULL
    if(length(res <- grep(".res", fefiles, value = TRUE))) {
      if(length(res2 <- res[!grepl("_df.res", res)]))
        for(tf in res2) {
          FixedEffects <- rbind(FixedEffects, df2m2(read.table(paste(dir, "/", tf, sep = ""), 
            header = TRUE)))
        }
      if(length(res2 <- res[grepl("_df.res", res)])) {
        nres2 <- length(res2) 
        fdf <- vector("list", length = nres2)
        for(k in 1L:nres2)
          fdf[[k]] <- read.table(paste(dir, "/", res2[k], sep = ""), header = TRUE)
      }
    }
    if(length(res <- grep("_sample.raw", fefiles, value = TRUE)))
      for(tf in res)
        fsample <- cbind(fsample, df2m(read.table(paste(dir, "/", tf, sep = ""), header = TRUE)))
    attr(FixedEffects, "sample") <- fsample
    attr(FixedEffects, "df") <- fdf
    rval$fixed.effects <- FixedEffects
  }

  ## create effect output
  if(!is.null(FixedEffects) && !is.null(data)) {
    rF <- rownames(FixedEffects)
    cF <- colnames(FixedEffects)
    FEattr <- attributes(FixedEffects)
    if(length(FixedEffects <- FixedEffects[rF != "(Intercept)",]))
      if(length(vars <- rF[rF %in% names(data)])) {
        j <- 0L
        if(is.vector(FixedEffects))
          FixedEffects <- matrix(FixedEffects, nrow = 1L)
        if("pstd" %in% cF) {
          id <- c(1L, 3L:ncol(FixedEffects))
          FixedEffects <- FixedEffects[,id]
          cF <- cF[id]
        }
        if(is.vector(FixedEffects))
          FixedEffects <- matrix(FixedEffects, nrow = 1L)
        colnames(FixedEffects) <- cF
        FEattrn <- names(FEattr)
        for(i in 1L:length(FEattrn))
          if(FEattrn[i] != "dim" && FEattrn[i] != "dimnames")
            attr(FixedEffects, FEattrn[i]) <- FEattr[[i]]
        for(tv in vars) {
          j <- j + 1L
          x <- unique(as.vector(unlist(data[tv])))
          vc <- matrix(FixedEffects[rownames(FixedEffects) == tv,], nrow = 1L)
          x <- cbind(x, x%*%vc)    
          x <- x[order(x[,1L]),]
          if(!is.matrix(x))
            x <- matrix(x, nrow = 1L)
          colnames(x) <- c(tv, cF)
          rownames(x) <- 1L:nrow(x)
          attr(x,"specs") <- list(dim = 1L, term = tv, label = tv)
          attr(x, "partial.resids") <- blow.up.resid(data, x, tv, 
            response, eta, 1L, "linear.bayesx")
          if(!is.null(attr(FixedEffects, "sample")))
            attr(x,"sample") <- attr(FixedEffects, "sample")[,j]
          class(x) <- "linear.bayesx"
          eval(parse(text = paste("rval$effects$\'", tv, "\' <- x", sep = "")))
        }
      }
  }

  return(rval)
}


unique.id <- function(x)
{
  rval <- .Call("unique_id",
    as.numeric(x),
    as.numeric(unique(x)))

  return(rval)
}


get.scale <- function(files, dir)
{
  var <- NULL
  if(any(grep("scale.res",files))) {
    sample <- NULL
    if(length(var <- grep("scale", files, value = TRUE))) {
      sc <- grep("sample", var)
      if(any(sc)) {
        sample <- grep("sample", var, value = TRUE)
        sample <- read.table(paste(dir, "/", sample, sep = ""), header = TRUE)
        sample$intnr <- NULL
        sample <- as.numeric(as.matrix(sample))
        id <- 1L:length(var)
        var <- var[id != sc]
      } else var <- var[1L]
    }
    var <- paste(dir, "/", var, sep = "")
    var <- read.table(var, header = TRUE)
    sn <- colnames(var)
    var <- as.matrix(var)
    if(length(sn) < 2L) {
      if(sn=="scale")
        sn <- "Scale"
      else
        sn <- "Sigma2"
      rownames(var) <- colnames(var) <- sn
    } else rownames(var) <- "Sigma2"
    attr(var,"sample") <- sample
  }

  return(var)
}


## reformate output
bayesx.reformate <- function(x)
{
  if(!is.null(x$fixed.effects) && length(x$fixed.effects) > 0L) {
    e <- x$fixed.effects
    type <- if(colnames(e)[1] == "pmode") "REML" else "MCMC"
    if(type == "REML") {
      tvalues <- e[,1L]/e[,4L]
      if(!is.null(x$model.fit$N)) {
        pvalues <- 2 * pt(-abs(tvalues), df = x$model.fit$N - 1)
        e <- cbind(e[,1L], e[,4L], tvalues, pvalues)
        if(!is.matrix(e))
          e <- matrix(e, nrow = 1L)
        colnames(e) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
      } else {
        e <- cbind(e[,1L], e[,4L], tvalues)
        if(!is.matrix(e))
          e <- matrix(e, nrow = 1L)
        colnames(e) <- c("Estimate", "Std. Error", "t value")
      }
    }
    if(type == "MCMC") {
      e <- e[, c(1L, 2L, 3L, 5L, 7L)]
      if(!is.matrix(e))
        e <- matrix(e, nrow = 1L)
      colnames(e) <- c("Mean", "Sd", "2.5%", "50%", "97.5%")
    }
    eattrn <- names(attributes(x$fixed.effects))
    for(i in 1L:length(eattrn))
      if(eattrn[i] != "dim" && eattrn[i] != "dimnames")
        attr(e, eattrn[i]) <- attr(x$fixed.effects, eattrn[i])
    rownames(e) <- rownames(x$fixed.effects)
    x$fixed.effects <- e
  }
  if(!is.null(x$effects) && length(x$effects) > 0L)
    for(i in 1L:length(x$effects)) {
      if(is.list(x$effects[[i]])) {
        if(length(x$effects[[i]]) > 0L)
          for(j in 1L:length(x$effects[[i]]))
            x$effects[[i]][[j]] <- chacol(x$effects[[i]][[j]])
      } else x$effects[[i]] <- chacol(x$effects[[i]])
    }
  if(!is.null(x$smooth.hyp) && length(x$smooth.hyp) > 0L)
    x$smooth.hyp <- recol(chacol(x$smooth.hyp))
  if(!is.null(x$random.hyp) && length(x$random.hyp) > 0L)
    x$random.hyp <- recol(chacol(x$random.hyp))
  if(!is.null(x$variance) && length(x$variance) > 0L)
    x$variance <- recol(chacol(x$variance))

  return(x)
}


chacol <- function(x)
{
  if(!is.null(cn <- colnames(x))) {
    for(i in 1L:length(cn)) {
      if(cn[i] == "pmean")
        cn[i] <- "Mean"
      if(cn[i] == "pqu2p5")
        cn[i] <- "2.5%"
      if(cn[i] == "pqu10")
        cn[i] <- "10%"
      if(cn[i] == "pmed")
        cn[i] <- "50%"
      if(cn[i] == "pqu50")
        cn[i] <- "50%"
      if(cn[i] == "pqu90")
        cn[i] <- "90%"
      if(cn[i] == "pqu97p5")
        cn[i] <- "97.5%"
      if(cn[i] == "pstddev")
        cn[i] <- "Sd"
      if(cn[i] == "pmin")
        cn[i] <- "Min"
      if(cn[i] == "pmax")
        cn[i] <- "Max"
      if(cn[i] == "pmode")
        cn[i] <- "Estimate"
      if(cn[i] == "ci95lower")
        cn[i] <- "2.5%"
      if(cn[i] == "ci80lower")
        cn[i] <- "10%"
      if(cn[i] == "std")
        cn[i] <- "std. error"
      if(cn[i] == "ci80upper")
        cn[i] <- "90%"
      if(cn[i] == "ci95upper")
        cn[i] <- "97.5%"
      if(cn[i] == "variance")
        cn[i] <- "Variance"
      if(cn[i] == "smoothpar")
        cn[i] <- "Smooth Par."
      if(cn[i] == "stopped")
        cn[i] <- "Stopped"
      if(cn[i] == "loglike")
        cn[i] <- "logLik"
      if(cn[i] == "aic")
        cn[i] <- "AIC"
      if(cn[i] == "bic")
        cn[i] <- "BIC"
      if(cn[i] == "gcv")
        cn[i] <- "GCV"
      if(cn[i] == "dic")
        cn[i] <- "DIC"
    }
    colnames(x) <- cn
  }

  return(x)
}


recol <- function(x)
{
  xattr <- attributes(x)
  xattrn <- names(xattr)
  rn <- rownames(x)
  ok <- FALSE
  if(any(id <- grepl("10%", colnames(x)))) {
    ok <- TRUE
    cn <- colnames(x)
    cn <- cn[!id]
    x <- x[,!id]
    if(!is.matrix(x)) {
      x <- matrix(x, nrow = 1L)
      rownames(x) <- rn
      colnames(x) <- cn
    }
  }
  if(any(id <- grepl("90%", colnames(x)))) {
    ok <- TRUE
    cn <- colnames(x)
    cn <- cn[!id]
    x <- x[,!id]
    if(!is.matrix(x)) {
      x <- matrix(x, nrow = 1L)
      rownames(x) <- rn
      colnames(x) <- cn
    }
  }
  if(ok)
    for(i in 1L:length(xattrn))
      if(xattrn[i] != "dim" && xattrn[i] != "dimnames")
        attr(x, xattrn[i]) <- xattr[[i]]

  return(x)
}
