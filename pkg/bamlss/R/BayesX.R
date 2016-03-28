######################
## BayesX interface ##
######################
BayesX.control <- function(n.iter = 1200, thin = 1, burnin = 200,
  seed = NULL, predict = "light", model.name = "bamlss", data.name = "d",
  prg.name = NULL, dir = NULL, verbose = FALSE, cores = NULL, ...)
{
  if(is.null(seed))
    seed <- '##seed##'
  stopifnot(burnin < n.iter)
  if(is.null(model.name))
    model.name <- 'bamlss'
  if(is.null(data.name))
    data.name <- 'd'
  if(is.null(prg.name))
    prg.name <- paste(model.name, 'prg', sep = '.')
  if(!grepl(".prg", prg.name))
    prg.name <- paste(prg.name, "prg", sep = ".")
  if(is.null(dir)) {
    dir.create(dir <- tempfile())
    attr(dir, "unlink") <- TRUE
  } else dir <- path.expand(dir)
  if(!file.exists(dir)) dir.create(dir)
  if(is.null(cores)) cores <- 1

  cvals <- list(
    "prg" = list(
      "iterations" = n.iter, "burnin" = burnin, "step" = thin,
      "setseed" = seed, "predict" = predict
    ),
    "setup" = list(
      "main" = c(rep(FALSE, 3), rep(TRUE, 2)), "model.name" = model.name, "data.name" = data.name,
      "prg.name" = prg.name, "dir" = dir, "verbose" = verbose, "cores" = cores
    )
  )

  cvals
}


BayesX <- function(x, y, family, start = NULL, weights = NULL, offset = NULL,
  data = NULL, control = BayesX.control(...), ...)
{
  if(is.null(family$bayesx))
    stop("BayesX specifications missing in family object, cannot set up model!")

  fbx <- family$bayesx

  model.name <- control$setup$model.name
  data.name <- rmf(control$setup$data.name)
  prg.name <- control$setup$prg.name
  dir <- control$setup$dir
  cores <- control$setup$cores
  if(is.null(cores)) cores <- 1

  if(!file.exists(dir)) {
    dir.create(dir)
    on.exit(unlink(dir))
  }

  if(is.null(data)) {
    data <- try(get('bf', envir = parent.frame())$model.frame, silent = TRUE)
    if(inherits(data, "try-error"))
      stop("cannot find the model.frame for creating BayesX model term objects!")
  }
  if(is.null(data))
    stop("no data available for creating BayesX model term objects!")

  yname <- colnames(y)[1]

  single_eqn <- function(x, y, id) {
    rhs <- dfiles <- prgex <- sdata <- NULL

    if(!is.null(x$model.matrix)) {
      cn <- rmf(colnames(x$model.matrix))
      colnames(x$model.matrix) <- cn
      cn <- cn[cn != "Intercept"]
      if(length(cn)) {
        rhs <- c(rhs, cn)
        sdata <- as.data.frame(x$model.matrix[, cn, drop = FALSE])
      }
      if("Intercept" %in% colnames(x$model.matrix)) {
        sdata <- cbind("Intercept" = rep(1, nrow(y)), sdata)
        rhs <- c("const", rhs)
      }
    }

    if(!is.null(sdata))
      sdata <- as.data.frame(sdata)

    if(!is.null(x$smooth.construct)) {
      for(j in names(x$smooth.construct)) {
        if(is.null(x$smooth.construct[[j]]$sx.construct))
          class(x$smooth.construct[[j]]) <- "userdefined.smooth.spec"
        sxc <- sx.construct(x$smooth.construct[[j]], data, id = c(id, j), dir = dir, mcmcreg = TRUE)
        if(!is.null(attr(sxc, "write")))
          prgex <- c(prgex, attr(sxc, "write")(dir))
        rhs <- c(rhs, sxc)
        tl <- x$smooth.construct[[j]]$term
        if(inherits(x$smooth.construct[[j]], "userdefined.smooth.spec"))
          tl <- paste(tl, collapse = "")
        if(!is.null(sdata)) {
          if(!all(tl %in% colnames(sdata))) {
            if(inherits(x$smooth.construct[[j]], "userdefined.smooth.spec")) {
              sdata[[tl]] <- runif(nrow(data))
            } else {
              sdata <- data[, tl, drop = FALSE]
            }
          }
        } else {
          if(inherits(x$smooth.construct[[j]], "userdefined.smooth.spec")) {
            sdata <- data[, 1, drop = FALSE]
            colnames(sdata) <- tl
            sdata[[1]] <- runif(nrow(data))
          } else {
            sdata <- data[, tl, drop = FALSE]
          }
        }
      }
    }

    rn <- response.name(as.formula(x$formula), hierarchical = FALSE)
    if(rn %in% family$names)
      rn <- NA
    if(is.na(rn))
      rn <- yname
    eqn <- paste(rn, "=", paste(rhs, collapse = " + "))
    rval <- list("eqn" = eqn, "prgex" = prgex)

    if(!is.null(sdata)) {
      if(nrow(sdata) == nrow(y))
        sdata <- cbind(sdata, y)
      rval$dname <- paste(paste(id, collapse = "_"), data.name, sep = "_")
      write.table(sdata, file = file.path(dir, paste(rval$dname, ".raw", sep = "")),
        quote = FALSE, row.names = FALSE)
      rval$prgex <- c(
        paste("dataset", rval$dname),
        paste(rval$dname, ".infile using ",
          file.path(dir, paste(rval$dname, ".raw", sep = "")), sep = ""),
        rval$prgex
      )
    }

    rval
  }

  eqn <- list()
  prgex <- NULL
  n <- 1
  for(i in names(x)) {
    if(!all(c("fake.formula", "formula") %in% names(x[[i]]))) {
      stop("hierarchical models not supported yet!")
      eqn[[i]] <- list()
      k <- 1
      for(j in names(x[[i]])) {
        msp <- single_eqn(x[[i]][[j]], y, id = c(i, j))
        teqn <- paste(model.name, ".hregress ", msp$eqn,
          ", family=", if(k < 2) fbx[[i]][1] else "gaussian_re",
          " equationtype=", fbx[[i]][2],
          if(n == length(x) & k < 2) {
            paste(" ", paste(names(control$prg), "=", control$prg, sep = "", collapse = " "))
          } else NULL,
          if(!is.null(msp$dname)) paste(" using", msp$dname) else NULL, sep = "")
        eqn[[i]][[j]] <- teqn
        prgex <- c(prgex, msp$prgex)
        k <- k + 1
      }
      eqn[[i]] <- rev(eqn[[i]])
    } else {
      msp <- single_eqn(x[[i]], y, id = i)
      teqn <- paste(model.name, ".hregress ", msp$eqn, ", family=", fbx[[i]][1],
        " equationtype=", fbx[[i]][2],
        if(n == length(x)) {
          paste(" ", paste(names(control$prg), "=", control$prg, sep = "", collapse = " "))
        } else NULL,
        if(!is.null(msp$dname)) paste(" using", msp$dname) else NULL, sep = "")
      eqn[[i]] <- teqn
      prgex <- c(prgex, msp$prgex)
    }
    n <- n + 1
  }

  prg <- c(prgex, "", paste("mcmcreg", model.name), "")
  for(i in unlist(rev(eqn)))
    prg <- c(prg, i, "")

  prg <- c(prg, paste(model.name, "getsample", sep = "."))
  prg <- c(
    paste('%% BayesX program created by bamlss: ', as.character(Sys.time()), sep = ''),
    paste('%% usefile ', file.path(dir, prg.name), sep = ''), "",
    prg
  )

  if(any(grepl("##seed##", prg, fixed = TRUE)))
    prg <- gsub("##seed##", round(runif(1L) * .Machine$integer.max), prg, fixed = TRUE)

#  cores <- control$setup$cores

#  if(cores > 1) {
#    k <- 1
#    while(file.exists(cdir <- file.path(dir, paste("Core", k, sep = "")))) {
#      k <- k + 1
#    }
#    dir.create(cdir)
#    prg <- gsub('##outfile##', prg <- file.path(cdir, model.name), prg, fixed = TRUE)
#    prg.name <- paste(gsub(".prg", "", prg.name, fixed = TRUE), k, ".prg", sep = "")
#    if(length(i <- grep('usefile', prg)))
#      prg[i] <- paste('%% usefile', file.path(dir, prg.name))
#    prg <- gsub(file.path(dir, paste(model.name, 'log', sep = '.')),
#      file.path(cdir, paste(model.name, 'log', sep = '.')), prg, fixed = TRUE)
#    prg <- gsub(paste(model.name, ".", sep = ""), paste(model.name, k, ".", sep = ""),
#      prg, fixed = TRUE)
#    prg[grep("mcmcreg", prg)] <- paste("mcmcreg ", model.name, k, sep = "")

#    dir <- cdir
#  }

  prgf <- file.path(dir, prg.name)
  writeLines(prg, prgf)

  require("BayesXsrc")

  warn <- getOption("warn")
  options(warn = -1)
  ok <- run.bayesx(prg = prgf, verbose = control$setup$verbose)
  options("warn" = warn)
  if(length(i <- grep("error", ok$log, ignore.case = TRUE))) {
    errl <- gsub("^ +", "", ok$log[i])
    errl <- gsub(" +$", "", errl)
    errl <- encodeString(errl, width = NA, justify = "left")
    errl <- paste(" *", errl)
    errm <- paste("an error occurred running the BayesX binary! The following messages are returned:\n",
      paste(errl, collapse = "\n", sep = ""), sep = "")
    warning(errm, call. = FALSE)
  }
  if(length(i <- grep("-nan", ok$log, ignore.case = TRUE))){
    warning("the BayesX engine returned NA samples, please check your model specification! In some cases it can be helpful to center continuous covariates!", call. = FALSE)
  }

  sfiles <- grep("_sample.raw", dir(file.path(dir, "output")), fixed = TRUE, value = TRUE)

  samples <- NULL
  for(i in names(x)) {
    if(!all(c("fake.formula", "formula") %in% names(x[[i]]))) {
      stop("hierarchical models not supported yet!")
    } else {
      if(!is.null(x[[i]]$model.matrix)) {
        sf <- grepl(paste("_", i, "_", sep = ""), sfiles, fixed = TRUE) & grepl(paste("LinearEffects", sep = ""), sfiles, fixed = TRUE)
        sf <- sfiles[sf]
        samps <- as.matrix(read.table(file.path(dir, "output", sf), header = TRUE)[, -1, drop = FALSE])
        colnames(samps) <- paste(i, ".p.", colnames(x[[i]]$model.matrix), sep = "")
        samples <- cbind(samples, samps)
      }
      if(!is.null(x[[i]]$smooth.construct)) {
        for(j in seq_along(x[[i]]$smooth.construct)) {
          tl <- x[[i]]$smooth.construct[[j]]$term
          if(is.null(x[[i]]$smooth.construct[[j]]$sx.construct))
            tl <- paste(tl, collapse = "")
          term <- paste("_", paste(tl, collapse = "_"), "_", sep = "")
          sf <- grepl(paste("_", i, "_", sep = ""), sfiles, fixed = TRUE) & grepl(term, sfiles, fixed = TRUE) & !grepl("_variance_", sfiles, fixed = TRUE)
          sf <- sfiles[sf]
          samps <- as.matrix(read.table(file.path(dir, "output", sf), header = TRUE)[, -1, drop = FALSE])
          cn <- colnames(x[[i]]$smooth.construct[[j]]$X)
          if(is.null(cn))
            cn <- paste("b", 1:ncol(x[[i]]$smooth.construct[[j]]$X), sep = "")
          colnames(samps) <- paste(i, ".s.", x[[i]]$smooth.construct[[j]]$label, ".", cn, sep = "")
          sf <- grepl(paste("_", i, "_", sep = ""), sfiles, fixed = TRUE) & grepl(term, sfiles, fixed = TRUE) & grepl("_variance_", sfiles, fixed = TRUE)
          sf <- sfiles[sf]
          if(length(sf)) {
            vsamps <- as.matrix(read.table(file.path(dir, "output", sf), header = TRUE)[, -1, drop = FALSE])
            colnames(vsamps) <- paste(i, ".s.", x[[i]]$smooth.construct[[j]]$label, ".", paste("tau2", 1:ncol(vsamps), sep = ""), sep = "")
            samps <- cbind(samps, vsamps)
          }
          samples <- cbind(samples, samps)
        }
      }
    }
  }

  if(FALSE) {
    sf <- grep("_DIC", dir(file.path(dir, "output")), fixed = TRUE, value = TRUE)
    dic <- read.table(file.path(dir, "output", sf), header = TRUE)[, -1, drop = FALSE]
    samples <- cbind(samples, "DIC" = dic$dic, "pd" = dic$pd)
  }

  as.mcmc(samples)
}


########################################
## (2) BayesX model term construction ##
########################################
sx <- function(x, z = NULL, bs = "ps", by = NA, ...)
{
  by <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
  call <- match.call()

  available.terms <- c(
    "rw1", "rw2",
    "season",
    "ps", "psplinerw1", "psplinerw2", "pspline",
    "te", "pspline2dimrw2", "te1", "pspline2dimrw1",
    "kr", "kriging",
    "gk", "geokriging",
    "gs", "geospline",
    "mrf", "spatial",
    "bl", "baseline",
    "factor",    
    "ridge", "lasso", "nigmix",
    "re", "ra", "random",
    "cs", "catspecific",
    "offset",
    "generic",
    "rps", "hrandom_pspline"
  )
  if(!bs %in% available.terms) stop(paste("basis type", sQuote(bs), "not supported by BayesX"))

  if(bs %in% c("rsps", "hrandom_pspline")) {
    bs <- "rsps"
    x <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
    rcall <- paste("r(x = ", by, ", bs = ", sQuote(bs), ", by = ", x, ", ...)", sep = "")
    rval <- eval(parse(text = rcall))
  } else {
    if(length(grep("~", term <- deparse(call$x))) && bs %in% c("re", "ra", "random")) {
      x <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
      rcall <- paste("r(x = ", x, ", by = ", by, ", ...)", sep = "")
      rval <- eval(parse(text = rcall))
    } else {
      k <- -1
      m <- NA
      xt <- list(...)
      if(is.null(xt$lambda))
        xt$lambda <- 100
      if("m" %in% names(xt))
        stop("argument m is not allowed, please see function s() using this specification!")
      if("k" %in% names(xt))
        stop("argument k is not allowed, please see function s() using this specification!")
      if(!is.null(xt$xt))
        xt <- xt$xt
      warn <- getOption("warn")
      options("warn" = -1)
      if(by != "NA" && is.vector(by) && length(by) < 2L && !is.na(as.numeric(by))) {
        xt["b"] <- by
        by <- "NA"
      }
      options("warn" = warn)
      if(bs %in% c("pspline2dimrw1", "pspline2dimrw2", "te",
        "gs", "geospline", "kr", "gk", "kriging", "geokriging")) {
        if(by != "NA")
          stop(paste("by variables are not allowed for smooths of type bs = '", bs, "'!", sep = ""))
      }
      if(bs == "te1") bs <- "pspline2dimrw1"
      if(!is.null(xt$knots))
        xt$nrknots <- xt$knots
      if(bs %in% c("ps", "te", "psplinerw1", "psplinerw2", "pspline",
        "pspline2dimrw2", "pspline2dimrw1", "gs", "geospline")) {
        if(!is.null(xt$degree))
          m <- xt$degree
        if(!is.null(xt$order)) {
          if(is.na(m))
            m <- c(3L, xt$order)
          else
            m <- c(m[1L], xt$order)
        }
        if(length(m) < 2L && (bs %in% c("gs", "geospline")))
          m <- c(3L, 1L)
        if(length(m) < 2L && is.na(m))
          m <- c(3L, 2L)   
        if(is.null(xt$order) && length(m) < 2L)
          m <- c(m, 2L)
        if(is.null(xt$order) && length(m) < 2L)
          m <- c(m, 1L)
        m[1L] <- m[1L] - 1L
        if(!is.null(xt$nrknots))
          k <- xt$nrknots + m[1L]
        else {
          if(bs %in% c("ps", "psplinerw1", "psplinerw2", "pspline"))
            k <- 20L + m[1L]
          else
            k <- 10L + m[1L]
        }
      }
      if(bs %in% c("kr", "gk", "kriging", "geokriging")) {
        m <- c(1L, 1L)
        if(!is.null(xt$nrknots)) {
          k <- xt$nrknots
        } else k <- -1L
      }
      if(!is.null(xt$map))
        xt$map.name <- rmf(as.character(call$map))
      xt[c("degree", "order", "knots", "nrknots")] <- NULL
      if(!length(xt))
        xt <- NULL
      if(!is.null(call$z)) 
        term <- c(term, deparse(call$z))
      if(bs != "te") {
        rval <- s(x, z, k = k, bs = bs, m = m, xt = xt)
      } else {
        rval <- te(x, z, k = k, bs = "ps", m = m, xt = xt, mp = FALSE)
        rval$margin[[1]]$term <- term[1]
        rval$margin[[2]]$term <- term[1]
      }
      rval$by <- by
      rval$term <- term
      rval$dim <- length(term)
      rval$label <- paste("sx(", paste(term, collapse = ",", sep = ""),
        if(by != "NA") paste(",by=", by, sep = "") else NULL, ")", sep = "")
    }
  }

  rval$special <- TRUE
  rval$sx.construct <- TRUE
  class(rval) <- c(class(rval), "no.mgcv")

  return(rval)
}

sx.construct <- function(object, data, ...)
{
  UseMethod("sx.construct")
}

sx.construct.default <- function(object, data, ...) 
{
  cl <- grep(".smooth.spec", class(object), value = TRUE, fixed = TRUE)
  bs <- gsub(".smooth.spec", "", cl, fixed = TRUE)
  info <- if(cl == paste(bs, "smooth.spec", sep = ".")) {
    paste("with basis", sQuote(bs))
  } else {
    paste("of class", sQuote(cl))
  }
  stop(paste("BayesX does not support smooth terms ", info,
    ", it is recommended to use sx() for specifying smooth terms",
    sep = ""))
}

do.xt <- function(term, object, not = NULL, noco = FALSE)
{
  if(!is.null(object$xt)) {
    names.xt <- names(object$xt)
    if(is.null(not))
      not <- "not"
    count <- 1
    co <- ","
    if(noco)
      co <- NULL
    for(name in names.xt) {
      if(count > 1)
        co <- ","
      if(!name %in% not) {
        if(name %in% c("full", "catspecific", "center", "derivative", "nofixed") || 
          is.logical(object$xt[[name]])) {
          if(is.logical(object$xt[[name]])) {
            if(object$xt[[name]])
              term <- paste(term, co, name, sep = "")
          } else term <- paste(term, co, name, sep = "")
          count <- count + 1
        } else {
          term <- paste(term, co, name, "=", object$xt[name], sep = "")
          count <- count + 1
        }
      }
    }
  }

  return(term)
}

sx.construct.userdefined.smooth.spec <- function(object, data, id, dir, ...)
{
  id <- paste(rmf(id), collapse = "_")
  term <- paste(object$term, collapse = "")
  Sn <- paste(id, "S", sep = "_")
  Xn <- paste(id, "X", sep = "_")
  if(length(object$S) > 1) {
    object$S <- list(do.call("+", object$S))
    object$rank <- qr(object$S[[1]])$rank
  }
  if(is.null(object$rank))
    object$rank <- qr(object$S[[1]])$rank
  if(is.null(object$xt$nocenter))
    object$xt$nocenter <- TRUE
  ##if(is.null(object$xt$centermethod))
  ##  object$xt$centermethod <- "meanfd"
  term <- paste(term, "(userdefined,penmatdata=", Sn, ",designmatdata=", Xn,
    ",rankK=", object$rank, sep = "")
  term <- paste(do.xt(term, object, c("center", "before")), ")", sep = "")

  write <- function(dir) {
    write.table(object$S[[1]], file = file.path(dir, paste(Sn, ".raw", sep = "")),
      quote = FALSE, row.names = FALSE)
    write.table(object$X, file = file.path(dir, paste(Xn, ".raw", sep = "")),
      quote = FALSE, row.names = FALSE)
    c(
      paste("dataset", Sn),
      paste(Sn, ".infile using ", file.path(dir, paste(Sn, ".raw", sep = "")), sep = ""),
      paste("dataset", Xn),
      paste(Xn, ".infile using ", file.path(dir, paste(Xn, ".raw", sep = "")), sep = "")
    )
  }

  attr(term, "write") <- write

  term
}

sx.construct.pspline.smooth <- sx.construct.ps.smooth.spec <- sx.construct.psplinerw1.smooth.spec <-
sx.construct.psplinerw2.smooth.spec <- sx.construct.pspline.smooth.spec <-
function(object, data, mcmcreg = FALSE, ...)
{
  if(length(object$p.order) == 1L) 
    m <- rep(object$p.order, 2L)
  else 
    m <- object$p.order
  m[is.na(m)] <- 2L
  object$p.order <- m
  object$p.order[1L] <- object$p.order[1L] + 1L
  if(inherits(object, "psplinerw1.smooth.spec"))
    object$p.order[2L] <- 1L
  if(inherits(object, "psplinerw2.smooth.spec"))
    object$p.order[2L] <- 2L
  if(object$bs.dim < 0L)
    object$bs.dim <- 10L
  if(length(object$p.order) > 1L) {
    if(object$p.order[2L] > 2L) {
      warning("order of the difference penalty not supported by BayesX, set to 2!")
      object$p.order <- c(object$p.order[1L], 2L)
    }
  }
  nrknots <- object$bs.dim - object$p.order[1L] + 1L
  if(nrknots < 5L) {
    warning("number of inner knots smaller than 5 not supported by BayesX, set to 5!",
      call. = FALSE)
    nrknots <- 5L
  }
  termo <- object$term
  term <- if(mcmcreg) {
    paste(termo, "(pspline,nrknots=",
      nrknots, ",degree=", object$p.order[1L], ",difforder=", object$p.order[2L], sep = "")
  } else {
    paste(termo, "(psplinerw", object$p.order[2L], ",nrknots=",
      nrknots, ",degree=", object$p.order[1L], sep = "")
  }
  object$xt[c("knots", "nrknots", "degree", "difforder")] <- NULL
  term <- paste(do.xt(term, object, c("center", "before")), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)
  
  return(term)
}

sx.construct.tensor.smooth <- sx.construct.tensor.smooth.spec <- sx.construct.t2.smooth.spec <-
function(object, dir, prg, data, mcmcreg = FALSE, ...)
{
  by <- object$term[1L]
  term <- object$term[2L]
  object <- object$margin[[1L]]
  object$term <- term
  object$by <- by
  term <- sx.construct(object, dir, prg, data, mcmcreg = mcmcreg, ...)
  if(!mcmcreg) {
    term <- gsub("psplinerw2", "pspline2dimrw2", term)
    term <- gsub("psplinerw1", "pspline2dimrw1", term)
  }
  return(term)
}

sx.construct.ra.smooth.spec <- sx.construct.re.smooth.spec <-
sx.construct.random.smooth.spec <- function(object, data, ...)
{
  term <- object$term
  if(is.null(object$ins))
    term <- paste(term, "(random", sep = "")
  else
    term <- paste(term, "(hrandom", sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}

sx.construct.rps.smooth.spec <- function(object, data, ...)
{
  term <- paste(object$term, "(hrandom_pspline,centermethod=meansum2", sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  term <- paste(object$by, term , sep = "*")

  return(term)
}

sx.construct.kr.smooth.spec <- sx.construct.kriging.smooth.spec <- function(object, data, ...)
{
  termo <- object$term
  if(length(termo) < 2L)
    stop("kriging method needs two terms!")
  if(object$bs.dim < 0L)
    object$bs.dim <- 10L
	nrknots <- object$bs.dim
  xt <- object$xt
  if(is.null(xt$full))
    term <- paste(termo[1L], "*", termo[2L], "(kriging,nrknots=", nrknots, sep = "")
  else {
    term <- paste(termo[1L], "*", termo[2L], "(kriging,full", sep = "")    
    object$xt$full <- NULL
  }
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}

construct.shrw <- function(object, data, what)
{
  term <- object$term
  term <- paste(term, "(", what, sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}

sx.construct.offset.smooth.spec <- function(object, data, ...)
{
  return(construct.shrw(object, data, "offset"))
}

sx.construct.mrf.smooth.spec <- sx.construct.spatial.smooth.spec <- function(object, data, ...)
{
  if(is.null(object$xt))
    stop("need to supply a map object in argument xt!")  
  map.name <- help.map.name(deparse(substitute(object, env = .GlobalEnv), 
    backtick = TRUE, width.cutoff = 500L))
  if(!is.null(object$xt$map.name))
    map.name <- object$xt$map.name
  if(!is.list(object$xt))
    object$xt <- list(object$xt)
  map.name <- rmf(gsub("\\s", "", paste(map.name, sep = "", collapse = "")))

  map <- object$xt$map
  if(is.null(map)) {
    if(!is.null(object$xt$polys))
      map <- object$xt$polys
    if(!is.null(object$xt$penalty))
      map <- object$xt$penalty
  }
  if(is.null(map))
    map <- object$xt$gra
  if(is.null(map)) {
    if(!is.list(object$xt[[1L]])) {
      if(inherits(object$xt[[1L]], "gra"))
        map <- object$xt[[1L]]
      else
        map <- object$xt
    } else map <- object$xt[[1L]]
    if(is.null(map)) {
      map <- object$xt
      if(is(map, "SpatialPolygonsDataFrame"))
        map <- SPDF2bnd(map)
      if(is.null(map) || (!is.list(map) && !inherits(map, "bnd") || !inherits(map, "gra")))
        stop("need to supply a bnd or graph file object in argument xt!")
    }
  }
  if(is(map, "nb"))
    map <- nb2gra(map)
  if(inherits(map, "SpatialPolygons"))
    map <- sp2bnd(map)
  if(!inherits(map, "bnd") && !inherits(map, "gra")) {
    if(is.list(map))
      class(map) <- "bnd"
    else
      class(map) <- "gra"
  }

  write <- function(dir = NULL) {
    if(is.null(dir)) {
      dir.create(dir <- tempfile())
      on.exit(unlink(dir))
    } else dir <- path.expand(dir)
    if(!file.exists(dir)) dir.create(dir)

    counter <- NULL
    ok <- TRUE
    files <- list.files(dir)
    while(ok) {
      classm <- class(map)
      if(length(classm) > 1L)
        if("list" %in% classm)
          class(map) <- classm[classm != "list"]
      mapfile <- paste(map.name, counter, ".", class(map), sep = "")[1]
      if(any(grepl(mapfile, files))) {
        if(is.null(counter))
          counter <- 0L
        counter <- counter + 1L
      } else ok <- FALSE
    }
    mapfile <- file.path(dir, mapfile)

    prg <- paste("map", map.name)
    if(inherits(map, "bnd")) {
      if(!any(is.na(poly.names <- as.integer(names(map))))) {
        poly.names <- sort(poly.names)
        poly.names <- as.character(poly.names)
      } else poly.names <- sort(names(map))
      map <- map[poly.names]
      class(map) <- "bnd"
      write.bnd(map = map, file = mapfile, replace = TRUE)
      prg <- c(prg, paste(map.name, ".infile using ", mapfile, sep = ""))
    } else {
      if(!is.character(map)) {
        write.gra(map = map, file = mapfile, replace = TRUE)
        prg <- c(prg, paste(map.name, ".infile, graph using ", mapfile, sep = ""))
      } else {
        stopifnot(is.character(map))
        pos <- regexpr("\\.([[:alnum:]]+)$", map)
        fext <- ifelse(pos > -1L, substring(map, pos + 1L), "")
        if(fext == "gra")
          prg <- c(prg, paste(map.name, ".infile, graph using ", path.expand(map), sep = ""))
        else
          prg <- c(prg, paste(map.name, ".infile using ", path.expand(map), sep = ""))
      }
    }
    prg
  }

  term <- object$term
  term <- paste(term, "(spatial,map=", map.name, sep = "")
  term <- paste(do.xt(term, object, c("map", "polys", "penalty", "map.name")), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  attr(term, "write") <- write
  attr(term, "map.name") <- map.name

  return(term)
}


make_by <- function(term, object, data)
{
  if(!missing(data) && !is.character(data)) {
    by <- data[[object$by]]
    if(is.factor(by) && nlevels(by) > 1) {
      nocenter <- paste(c("lasso", "nigmix", "ridge", "ra"), ".smooth.spec", sep = "")
      term <- paste(paste(rmf(object$by), rmf(levels(by)), sep = ""), "*", term, sep = "")
      if((k <- length(term)) > 1)
        for(j in 1:k) 
          if(!grepl("center", term[j]) && is.null(object$xt$center))
            if(!class(object) %in% nocenter)
              term[j] <- gsub(")", ",center)", term[j])
    } else term <- paste(rmf(object$by), "*", term, sep = "")
  } else term <- paste(rmf(object$by), "*", term, sep = "")

  return(term)
}

rmf <- function(x) 
{
  for(i in 1L:length(x)) {
    for(char in c("+", "-", "*", ":", "^", "/", " ", "(", ")", "]", "[",
      ",", ".", "<", ">", "?", "!", "'", "#", "~", "`", ";", "=", "&", "$", "@")) {
      x[i] <- gsub(char, "", x[i], fixed = TRUE)
    }
  }

  return(x)
}

help.map.name <- function(x)
{
  if(is.null(x))
    return("")
  x <- splitme(x)
  if(resplit(x[1L:2L]) == "s(") {
    x <- resplit(c("sfoofun", x[2L:length(x)]))
    x <- eval(parse(text = x), envir = parent.frame())
    if(is.null(x))
      x <- "map"
  } else  x <- "map"

  return(x)
}

sfoofun <- function(x, xt = NULL, ...)
{
  if(is.null(x) || is.null(xt))
    return(NULL)
  x <- rval <- deparse(substitute(xt), backtick = TRUE, width.cutoff = 500L)
  if(is.list(xt) && length(xt)>1)
    for(i in 1L:length(xt))
      if(inherits(xt[[i]], "bnd") || inherits(xt[[i]], "gra") || inherits(xt[[i]], "list")) {
        rval <- strsplit(x, ",", " ")[[1L]]
        if(length(rval) > 1L)
          rval <- rval[i]
      }
  rval <- splitme(rval)
  if(length(rval) > 5L)
    if(resplit(rval[1L:5L]) == "list(")
      rval <- rval[6L:length(rval)]
  if(rval[length(rval)] == ")")
    rval <- rval[1L:(length(rval) - 1)]
  if(any(grepl("=", rval)))
    rval <- rval[(grep("=", rval) + 2L):length(rval)]
  rval <- resplit(rval)
   
  return(rval)
}

splitme <- function(x) {
  return(strsplit(x, "")[[1L]])
}

resplit <- function(x) {
  if(!is.null(x))
    x <- paste(x, sep = "", collapse = "")
	
  return(x)
}

