## prediction method based on refitting with weights
predict.bayesx <- function(object, newdata = NULL, model = NULL,
  type = c("response", "link", "terms", "model"), na.action = na.pass, ...)
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
    trms <- terms(object)
    response <- as.character(trms)[2]
    newdata[[response]] <- 0
    newdata <- model.frame.bayesx(trms, data = newdata, na.action = na.action, ...)
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
    nd[[response]] <- c(mf[[response]], rep(NA, length = nrow(newdata)))
    nd <- as.data.frame(nd)
    weights <- model.weights(mf)
    if(is.null(weights))
      weights <- rep(1, length = nrow(mf))
    weights <- c(weights, rep(0, length = nrow(newdata)))
    nm <- update(object, . ~ ., data = nd, weights = weights,
      seed = object$bayesx.setup$setseed, parse.model.frame = FALSE)
    return(nm)
  }
  return(fitted.bayesx(object, ...))
}

