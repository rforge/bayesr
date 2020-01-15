jm2_bamlss <- function(k = 1)
{
  if(k > 0) {
    mu <- paste("mu", 1:k, sep = "")
    sigma <- paste("sigma", 1:k, sep = "")
    alpha <- paste("alpha", 1:k, sep = "")
    links <- c(rep("identity", k), rep("log", k), rep("identity", k), "identity")
    names(links) <- c(mu, sigma, alpha, "psi")
  }

  links <- c(lambda = "log", gamma = "log", links)

  rval <- list(
    "family" = "jm",
    "names" = names(links),
    "links" = links,
    "transform" = function(x, obstime = NULL, id = NULL, kpsi = 5, ...) {
      jm2_transform(x = x$x, y = x$y, terms = x$terms, knots = x$knots,
        formula = x$formula, family = x$family, data = x$model.frame,
        obstime = obstime, id = id, kpsi = kpsi, ...)
    }
  )
  
  class(rval) <- "family.bamlss"
  rval
}


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
  object$neigen <- object$bs.dim
  for(j in 1:object$neigen) {
    object$sfun[[j]] <- splinefun(object$fpca$workGrid,
      object$fpca$phi[, j] * sqrt(object$fpca$lambda[j]))
    Z[[j]] <- object$sfun[[j]](data[[object$term[2]]])
  }
  Z <- do.call("cbind", Z)
  object$form <- as.formula(paste("~", paste(object$term[1], collapse = ":"), "-1"))
  Re <- model.matrix(object$form, data)
  object$X <- tensor.prod.model.matrix(list(Re, Z))
  object$S <- list(diag(ncol(object$X)))

  Cr <- rbind(t(colSums(Re)), matrix(0, nrow = ncol(Re) - 1, ncol = ncol(Re)))

  object$bs.dim <- ncol(object$X)
  object$S <- list(diag(object$bs.dim))
  object$rank <- object$bs.dim
  object$null.space.dim <- 0
  object$C <- matrix(0, 0, ncol(object$X))
  object$Cr <- Cr
  object$side.constrain <- FALSE
  object$plot.me <- TRUE
  object$te.ok <- 2

  class(object) <- c("fri.smooth", "random.effect")
  return(object)
}

Predict.matrix.fri.smooth <- function(object, data)
{
  Z <- list()
  for(j in 1:object$neigen) {
    Z[[j]] <- object$sfun[[j]](data[[object$term[2]]])
  }
  Z <- do.call("cbind", Z)
  Re <- model.matrix(object$form, data)
  X <- tensor.prod.model.matrix(list(Re, Z))
  return(X)
}

jm2_transform <- function(x, y, terms, knots,
  formula, family, data, obstime, id, kpsi, ...)
{
  mu <- grep("mu", names(formula), value = TRUE)
  ## if length(mu) == 1, use formula from mu, else use multivariate fpca.
  if(length(mu)) {
    ynames <- NULL
    for(j in mu) {
      ynames <- c(ynames, response.name(formula[[j]]$formula))
      b <- gam(formula[[j]]$formula, data = d)
plot(b)
    }
print(ynames)
    Y <- data[, ynames, drop = FALSE]
    vtime <- data[[obstime]]
    vid <- data[[id]]
print(head(Y))
print(str(formula))
  }
  stop()
}

if(FALSE) {
  x <- sort(runif(100, 0, 3))
  y <- sin(x)

  sc <- smooth.construct(s(x,bs="ps"),list(x=x),NULL)

  beta <- coef(lm(y ~ -1 + sc$X))

  fit <- drop(sc$X %*% beta)

  plot(x, y)
  lines(fit ~ x, col = 2)

  library("statmod")

  gq <- gauss.quad(15)
  
  a <- min(x)
  b <- max(x)
  x2 <- gq$nodes * (b - a) / 2 + (a + b) / 2

  B <- PredictMat(sc, data.frame(x = x2))

  fit2 <- drop(B %*% beta)

  plot(fit2 ~ x2, type = "l", xlim = range(x))

  int <- (b - a) / 2 * sum(fit2 * gq$weights)

  f <- splinefun(x, y)
  integrate(f, a, b)
}

if(FALSE) {
library(fdapace)

################################################################################
# Functions to help the process
################################################################################
# Function to transform the time
t_transform <- function(var, data) {
  
  # awkward code for rescaling the observed censoring/event time
  cutted <- cut(data[, var], 
                breaks = c(0, 7, 14, 21, 28, 90, 180, 365, 746),
                # 746 = max(vall$survtime, na.rm = TRUE) + 1), 
                include.lowest = FALSE, right = FALSE,
                labels = paste0(1:8))
  myfun <- function(category, x){
    switch(as.character(category), 
           "1" = x / (7 - 0),
           "2" = (x - 7) / (14 - 7) + 1,
           "3" = (x - 14) / (21 - 14) + 2,
           "4" = (x - 21) / (28 - 21) + 3,
           "5" = (x - 28) / (90 - 28) + 4,
           "6" = (x - 90) / (180 - 90) + 5,
           "7" = (x - 180) / (365 - 180) + 6,
           "8" = (x - 365) / (746 - 365) + 7,
           "NA" = NA)
  }
  unlist(mapply(myfun, category = cutted, x = data[, var]))
  
}
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------

################################################################################
# Working with the data
################################################################################

load(file = "../Data/vall.Rdata")

# Prepare the data
vall$idpseud <- as.factor(vall$idpseud)

# Transform the time
vall$surv_t <- t_transform(var = "survtime", data = vall)
vall$obs_t <- t_transform(var = "obstime", data = vall)

# Transform the markers
for(j in grep("marker", names(vall)))
  vall[[j]] <- log(vall[[j]] + 1)

d <- with(vall, data.frame("m1" = marker1, "m2" = marker2, "id" = idpseud,
  "otime" = obs_t, "stime" = surv_t, "event" = event))

d <- na.omit(d)

i <- table(d$id)
i <- names(i)[i > 2]
d <- subset(d, id %in% i)
d$id <- droplevels(d$id)

f <- list(
  Surv(stime, event) ~ s(stime),
  gamma ~ 1,
  m1 ~ s(otime),
  m2 ~ s(otime),
  sigma1 ~ 1,
  sigma2 ~ 1,
  alpha1 ~ 1,
  alpha2 ~ 1,
  psi ~ otime + id
)

b <- bamlss(f, data = d, family = jm2_bamlss(k = 2), sampler = FALSE,
  obstime = "otime", id = "id")

}
