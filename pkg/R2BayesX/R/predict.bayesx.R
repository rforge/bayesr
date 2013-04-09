## prediction method based on refitting with weights
predict.bayesx <- function(object, newdata = NULL, model = NULL,
  type = c("response", "link", "terms", "model"), na.action = na.pass,
  digits = 5, ...)
{
  type <- match.arg(type)
  if(!is.null(newdata)) {
    object <- get.model(object, model)
    if(length(object) > 1)
      stop("argument model is specified wrong, predictions are only possible one single model!")

    object <- object[[1]]
    mf <- model.frame(object)
    stopifnot(inherits(newdata, c("data.frame", "list", "matrix")))
    newdata <- as.data.frame(newdata)

    ff <- object$formula
    response <- all.vars(ff)[1]
    newdata[[response]] <- 0
    newdata <- model.frame.bayesx(ff, data = newdata, na.action = na.action, ...)
    newdata[[response]] <- NULL

    nam_nd <- names(newdata)
    nam_mf <- names(mf)
    nam_mf <- nam_mf[nam_mf != response]
    if(!all(nc <- nam_mf %in% nam_nd))
      stop(paste("variables", paste(nam_mf[!nc], collapse = ", "), "are missing in newdata!"))

    nd <- list()
    for(j in nam_mf) {
      nd[[j]] <- c(mf[[j]], newdata[[j]])
    }
    nd <- as.data.frame(nd)
    names(nd) <- nam_mf
    nd[[response]] <- c(mf[[response]], rep(0, length = nrow(newdata)))

    weights <- model.weights(mf)
    if(is.null(weights))
      weights <- rep(1, length = nrow(mf))
    i <- c(rep(FALSE, length(weights)), rep(TRUE, nrow(newdata)))
    nd$weights <- c(weights, rep(0, length = nrow(newdata)))

    pm <- update(object, . ~ ., data = nd, weights = weights,
      seed = object$bayesx.setup$setseed, parse.model.frame = FALSE)
    
    if(type == "model") return(pm)
    if(type == "response") {
      pr <- fitted(pm)
      if(!is.null(dim(pr)))
        pr <- pr[, 1]
      pr <- pr[i]
    }
    if(type == "link") {
      pr <- fitted(pm)
      if(!is.null(dim(pr))) {
        if(any(j <- grepl("mu", names(pr))))
          pr <- pr[, j, drop = FALSE]
        pr <- pr[i, ]
      } else pr <- pr[i]
    }
    if(type == "terms") {
      pr <- fitted(pm, term = names(object$effects))
      if(inherits(pr, "data.frame")) pr <- list(pr)
      labels <- NULL
      for(j in seq_along(pr)) {
        if(inherits(pr[[j]], "data.frame")) {
          nt <- names(pr[[j]])[1]
          i <- round(pr[[j]][[nt]], digits) %in% round(newdata[[nt]], digits)
          pr[[j]] <- pr[[j]][i, ]
          rownames(pr[[j]]) <- NULL
        }
        tl <- attr(pr[[j]], "specs")$label
        labels <- c(labels, if(is.null(tl)) "NA" else tl)
      }
      names(pr) <- labels
      if(length(pr) < 2) pr <- pr[[1]]
    }
    return(pr)
  }
  return(fitted.bayesx(object, ...))
}

