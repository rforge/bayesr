# Function to transform data set to FunData objects
marker_to_irregFunData <- function (marker, data, t, id = "idpseud", 
                                    fundata = FALSE) {
  
  # + marker +
  # char string giving the name of the marker to be extracted
  # + data +
  # dataframe with all the markers
  # + id +
  # char string giving the name of the id variable
  # + t +
  # char string giving the name of the time variable
  data <- as.data.frame(data)
  # Exclude observations with missing marker values or time from the FunData 
  # object
  data$full_info <- !(is.na(data[, marker]) | is.na(data[, t]))
  
  # Extract time arguments
  argvals <- lapply(unique(data[, id]), function (x) {
    data[data[, id] == x & data$full_info, t]
  })
  
  # Extract marker values
  X <- lapply(unique(data[, id]), function (x) {
    data[data[, id] == x & data$full_info, marker]
  })
  
  if(fundata) {
    require(funData)
    fun <- irregFunData(argvals = argvals, X = X)
  } else {
    fun <- list(argvals = argvals, X = X)
  }

  fun
}


## Smooth constructor.
smooth.construct.fri.smooth.spec <- function(object, data, knots, ...)
{
  stopifnot(requireNamespace("fdapace"))
  if(is.null(object$xt$m))
    stop("no markers specified!")
  if(is.null(dim(object$xt$m)))
    object$xt$m <- list("m1" = object$xt$m)
  nm <- names(object$xt$m)
  if(is.null(nm)) {
    nm <- paste0("m", 1:length(object$xt$m))
    names(object$xt$m) <- nm
  }
  for(j in nm)
    data[[j]] <- object$xt$m[[j]]
  m <- marker_to_irregFunData(marker = nm, data = data,
    id = object$term[1], t = object$term[2])
  object$fpca <- fdapace::FPCA(m$X, m$argvals)
  object$sfun <- list()
  if(object$bs.dim < 1)
    object$bs.dim <- ncol(object$fpca$phi)
  Z <- list()
  for(j in 1:object$bs.dim) {
    object$sfun[[j]] <- splinefun(object$fpca$workGrid, object$fpca$phi[, j])
    Z[[j]] <- object$sfun[[j]](data[[object$term[2]]])
  }
  Z <- do.call("cbind", Z)
  object$form <- as.formula(paste("~", paste(object$term[1], collapse = ":"), "-1"))
  Re <- model.matrix(object$form, data)
  object$X <- tensor.prod.model.matrix(list(Re, Z))
  object$S <- list(diag(ncol(object$X)))
  object$rank <- qr(object$S[[1]])$rank
  object$null.space.dim <- ncol(object$S[[1]])
  class(object) <- "fri.smooth"
  return(object)
}

Predict.matrix.fri.smooth <- function(object, data)
{
  Z <- list()
  for(j in 1:object$bs.dim) {
    Z[[j]] <- object$sfun[[j]](data[[object$term[2]]])
  }
  Z <- do.call("cbind", Z)
  Re <- model.matrix(object$form, data)
  X <- tensor.prod.model.matrix(list(Re, Z))
  return(X)
}

