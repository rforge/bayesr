find.fixed.effects <-
function(dir, files, data, response, eta, model.name, rval, minfo)
{
  fixed2table <- function(file, tf) {
    tb <- readLines(file)
    tb <- gsub("   ", ",", tb)
    con <- paste(tempdir(), "/temp.data.raw", sep = "")
    tb[1] <- gsub(" ", ",", tb[1])
    stb <- strsplit(tb, "")
    for(i in 1:length(stb)) {
      if(stb[[i]][k <-length(stb[[i]])] == ",") 
        stb[[i]] <- stb[[i]][1:(k - 1)]
      tb[i] <- paste(stb[[i]], sep = "", collapse = "")
    }
    tb <- gsub("const", "(Intercept)", tb, fixed = TRUE)
    writeLines(tb, con)
    x <- read.table(con, header = TRUE, sep = ",")
    if(!is.null(vn <- as.character(x$varname))) {
      if(any(id <- grepl(" (cat. ", vn, fixed = TRUE))) {
        tmp <- strsplit(vn[id], "cat. ")
        new <- NULL
        for(split in tmp) {
          one <- strsplit(split[1L], " ")[[1L]][1L]
          two <- splitme(split[2L])
          two <- resplit(two[1:(length(two) - 1)])
          if(!is.null(minfo)) {
            if(!is.null(minfo$YLevels) && !is.null(minfo$nYLevels)) {
              oL <- eval(parse(text = minfo$YLevels))
              nL <- eval(parse(text = minfo$nYLevels))
              two <- oL[nL == two]
            }
          }
          new <- c(new, paste(one, ":", two, sep = ""))
        }
        vn[id] <- new
        x$varname <- vn
      }
    }
    if(NROW(x) == 1L && length(split <- strsplit(tf, "_")[[1L]]) == 3L) {
      if(!is.null(x$varname)) {
        split <- gsub(".res", "", split[3L])
        if(!is.null(minfo)) {
          if(!is.null(minfo$YLevels) && !is.null(minfo$nYLevels)) {
            oL <- eval(parse(text = minfo$YLevels))
            nL <- eval(parse(text = minfo$nYLevels))
            split <- oL[nL == split]
          }
        }
        x$varname <- paste(x$varname, ":", split, sep = "", collapse = "")
      }
    }
    return(x)
  }
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
        for(tf in res2)
          FixedEffects <- rbind(FixedEffects, df2m2(fixed2table(paste(dir, "/", tf, sep = ""), tf)))
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
        if(!is.matrix(FixedEffects)) {
          FixedEffects <- matrix(FixedEffects, nrow = 1)
          colnames(FixedEffects) <- cF
          rownames(FixedEffects) <- rF[rF != "(Intercept)"]
        }
        if("pstd" %in% cF) {
          id <- c(1L, 3L:ncol(FixedEffects))
          FixedEffects <- FixedEffects[,id]
          cF <- cF[id]
        }
        if(!is.matrix(FixedEffects)) {
          FixedEffects <- matrix(FixedEffects, nrow = 1)
          colnames(FixedEffects) <- cF
          rownames(FixedEffects) <- rF[rF != "(Intercept)"]
        }
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
          if(!is.null(attr(x, "partial.resids"))) {
            attr(x, "partial.resids")[,2L] <- attr(x, "partial.resids")[,2L] - 
              mean(attr(x, "partial.resids")[,2L])
          }
          if(!is.null(attr(FixedEffects, "sample")))
            attr(x,"sample") <- attr(FixedEffects, "sample")[,j]
          class(x) <- c("linear.bayesx", "matrix")
          eval(parse(text = paste("rval$effects$\'", tv, "\' <- x", sep = "")))
        }
      }
  }

  return(rval)
}

