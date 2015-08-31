## From mgcv.
expand.t2.smooths <- function(sm) 
{
  m <- length(sm)
  not.needed <- TRUE
  for(i in 1:m) if (inherits(sm[[i]], "t2.smooth") && length(sm[[i]]$S) > 1) {
    not.needed <- FALSE
    break
  }
  if(not.needed) 
    return(NULL)
  smr <- list()
  k <- 0
  for(i in 1:m) {
    if(inherits(sm[[i]], "t2.smooth")) {
      smi <- split.t2.smooth(sm[[i]])
      comp.ind <- (k + 1):(k + length(smi))
      for(j in 1:length(smi)) {
        k <- k + 1
        smr[[k]] <- smi[[j]]
        smr[[k]]$comp.ind <- comp.ind
      }
    } else {
      k <- k + 1
      smr[[k]] <- sm[[i]]
    }
  }
  smr
}

split.t2.smooth <- function (object) 
{
  if(!inherits(object, "t2.smooth")) 
    return(object)
  ind <- 1:ncol(object$S[[1]])
  ind.para <- object$first.para:object$last.para
  sm <- list()
  sm[[1]] <- object
  St <- object$S[[1]] * 0
  for(i in 1:length(object$S)) {
    indi <- ind[diag(object$S[[i]]) != 0]
    label <- paste(object$label, ".frag", i, sep = "")
    sm[[i]] <- list(S = list(object$S[[i]][indi, indi]), 
      first.para = min(ind.para[indi]), last.para = max(ind.para[indi]), 
      fx = object$fx[i], fixed = object$fx[i], sp = object$sp[i], 
      null.space.dim = 0, df = length(indi), rank = object$rank[i], 
      label = label, S.scale = object$S.scale[i])
    class(sm[[i]]) <- "t2.frag"
    St <- St + object$S[[i]]
  }
  i <- length(object$S) + 1
  indi <- ind[diag(St) == 0]
  if(length(indi)) {
    label <- paste(object$label, ".frag", i, sep = "")
    sm[[i]] <- list(S = NULL, first.para = min(ind.para[indi]), 
      last.para = max(ind.para[indi]), fx = TRUE, fixed = TRUE, 
      null.space.dim = 0, label = label, df = length(indi))
    class(sm[[i]]) <- "t2.frag"
  }
  sm
}


smooth2random <- function(object, vnames, type = 1) {
  UseMethod("smooth2random")
}

smooth2random.t2.smooth <- function (object, vnames, type = 1) 
{
    if (object$fixed) 
        return(list(fixed = TRUE, Xf = object$X))
    fixed <- rep(TRUE, ncol(object$X))
    random <- list()
    diagU <- rep(1, ncol(object$X))
    ind <- 1:ncol(object$X)
    pen.ind <- ind * 0
    n.para <- 0
    for (i in 1:length(object$S)) {
        indi <- ind[diag(object$S[[i]]) != 0]
        pen.ind[indi] <- i
        X <- object$X[, indi, drop = FALSE]
        D <- diag(object$S[[i]])[indi]
        diagU[indi] <- 1/sqrt(D)
        X <- X %*% diag(diagU[indi])
        fixed[indi] <- FALSE
        term.name <- new.name("Xr", vnames)
        group.name <- new.name("g", vnames)
        vnames <- c(vnames, term.name, group.name)
        if (type == 1) {
            form <- as.formula(paste("~", term.name, "-1", sep = ""), 
                env = .GlobalEnv)
            random[[i]] <- pdIdnot(form)
            names(random)[i] <- group.name
            attr(random[[i]], "group") <- factor(rep(1, nrow(X)))
            attr(random[[i]], "Xr.name") <- term.name
            attr(random[[i]], "Xr") <- X
        }
        else {
            random[[i]] <- X
            names(random)[i] <- term.name
            attr(random[[i]], "s.label") <- object$label
        }
        n.para <- n.para + ncol(X)
    }
    if (sum(fixed)) {
        Xf <- object$X[, fixed, drop = FALSE]
    }
    else Xf <- matrix(0, nrow(object$X), 0)
    list(rand = random, trans.D = diagU, Xf = Xf, fixed = FALSE, 
        rind = 1:n.para, rinc = rep(n.para, n.para), pen.ind = pen.ind)
}


smooth2random.mgcv.smooth <- function (object, vnames, type = 1) 
{
    if (object$fixed) 
        return(list(fixed = TRUE, Xf = object$X))
    if (length(object$S) > 1) 
        stop("Can not convert this smooth class to a random effect")
    ev <- eigen(object$S[[1]], symmetric = TRUE)
    null.rank <- object$df - object$rank
    p.rank <- object$rank
    if (p.rank > ncol(object$X)) 
        p.rank <- ncol(object$X)
    U <- ev$vectors
    D <- c(ev$values[1:p.rank], rep(1, null.rank))
    D <- 1/sqrt(D)
    UD <- t(t(U) * D)
    X <- object$X %*% UD
    if (p.rank < object$df) 
        Xf <- X[, (p.rank + 1):object$df, drop = FALSE]
    else Xf <- matrix(0, nrow(object$X), 0)
    term.name <- new.name("Xr", vnames)
    if (type == 1) {
        form <- as.formula(paste("~", term.name, "-1", sep = ""), 
            env = .GlobalEnv)
        random <- list(pdIdnot(form))
        group.name <- new.name("g", vnames)
        names(random) <- group.name
        attr(random[[1]], "group") <- factor(rep(1, nrow(X)))
        attr(random[[1]], "Xr.name") <- term.name
        attr(random[[1]], "Xr") <- X[, 1:p.rank, drop = FALSE]
    }
    else {
        random <- list(X[, 1:p.rank, drop = FALSE])
        names(random)[1] <- term.name
        attr(random[[1]], "s.label") <- object$label
    }
    rind <- 1:p.rank
    pen.ind <- rep(0, ncol(object$X))
    pen.ind[rind] <- 1
    rinc <- rep(p.rank, p.rank)
    list(rand = random, Xf = Xf, trans.U = U, trans.D = D, fixed = FALSE, 
        rind = rind, rinc = rinc, pen.ind = pen.ind)
}

smooth2random.tensor.smooth <- function (object, vnames, type = 1) 
{
    if (type == 2) 
        stop("te smooths not useable with gamm4: use t2 instead")
    if (sum(object$fx) == length(object$fx)) 
        return(list(fixed = TRUE, Xf = object$X))
    else if (sum(object$fx) != 0) 
        warning("gamm can not fix only some margins of tensor product.")
    sum.S <- object$S[[1]]/mean(abs(object$S[[1]]))
    if (length(object$S) > 1) 
        for (l in 2:length(object$S)) {
            sum.S <- sum.S + object$S[[l]]/mean(abs(object$S[[l]]))
        }
    null.rank <- object$null.space.dim
    ev <- eigen(sum.S, symmetric = TRUE)
    p.rank <- ncol(object$X) - null.rank
    if (p.rank > ncol(object$X)) 
        p.rank <- ncol(object$X)
    U <- ev$vectors
    D <- c(ev$values[1:p.rank], rep(1, null.rank))
    if (sum(D <= 0)) 
        stop("Tensor product penalty rank appears to be too low: please email Simon.Wood@R-project.org with details.")
    U <- U
    X <- object$X %*% U
    if (p.rank < ncol(X)) 
        Xf <- X[, (p.rank + 1):ncol(X), drop = FALSE]
    else Xf <- matrix(0, nrow(X), 0)
    for (l in 1:length(object$S)) {
        object$S[[l]] <- (t(U) %*% object$S[[l]] %*% U)[1:p.rank, 
            1:p.rank]
        object$S[[l]] <- (object$S[[l]] + t(object$S[[l]]))/2
    }
    term.name <- new.name("Xr", vnames)
    form <- as.formula(paste("~", term.name, "-1", sep = ""), 
        env = .GlobalEnv)
    attr(form, "S") <- object$S
    random <- list(pdTens(form))
    group.name <- new.name("g", vnames)
    names(random) <- group.name
    attr(random[[1]], "group") <- factor(rep(1, nrow(X)))
    attr(random[[1]], "Xr.name") <- term.name
    attr(random[[1]], "Xr") <- X[, 1:p.rank, drop = FALSE]
    rind <- 1:p.rank
    rinc <- rep(p.rank, p.rank)
    list(rand = random, Xf = Xf, trans.U = U, trans.D = rep(1, 
        ncol(U)), fixed = FALSE, rind = rind, rinc = rinc)
}


## From R2BayesX.
r <- function(x, h = NULL, by = NA, xt = NULL, 
  data = NULL, weights = NULL, subset = NULL, 
  offset = NULL, na.action = na.fail, contrasts = NULL, 
  control = bayesx.control(...), ...)
{
  term <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500L)
  call <- match.call()
  is.formula <- FALSE
  if(!any(grepl("~", term)) && is.null(h) && !is.null(data)) {
    term <- paste(term, "~", 1)
  } 
  if(any(grepl("~", term))) {
    tmp <- strsplit(term, "~")[[1L]]
    call$h <- h <- term
    term <- splitme(tmp[1L])
    term <- resplit(term[term != " "])
    call$x <- term
    is.formula <- TRUE
  }
  by.var <- if(!is.character(by)) deparse(substitute(by), backtick = TRUE, width.cutoff = 500L) else by
  ins <- formula <- NULL
  if(by.var == ".") 
    stop("by=. not allowed")
  if(term == ".") 
    stop("r(.) not yet supported.")
  label <- paste("r(", term)
  if(!is.null(h)) {
    ins <- list()
    mlabel <- paste(as.character(call$h), collapse = " ")
    split <- splitme(mlabel)
    if(split[1L] != "~" && !is.formula)
      mlabel <- resplit(c("~", split))
    if(!is.formula)
      formula <- as.formula(paste(term, mlabel))
    else
      formula <- as.formula(mlabel)
    if(length(grep("~", mlabel, fixed = TRUE)))
      label <- paste("r(", mlabel)
    else
      label <- paste(label, ",", mlabel, collapse="")
    mf <- terms.formula(formula, specials=c("sx", "s", "te", "r"))
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
  if(!is.null(control$bs) && control$bs == "rsps") {
    rval$control$bs <- NULL
    class(rval) <- "rsps.smooth.spec"
    if(is.null(rval$ins)) {
      rval$ins <- list()
      rval$formula <- as.formula(paste(rval$term, "~ -1"))
    }
  } else class(rval) <- "ra.smooth.spec"

  return(rval) 
}

