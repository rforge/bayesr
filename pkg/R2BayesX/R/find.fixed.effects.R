find.fixed.effects <-
function(dir, files, data, response, eta, model.name, rval)
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

