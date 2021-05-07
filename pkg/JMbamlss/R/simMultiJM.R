

# New Simulation Function For Multivariate JMs Based On FPCs --------------


# Adapt the structure given by simJM function in bamlss

#' @param nsub Number of subjects.
#' @param times Vector of time points.
#' @param probmiss Probability of missingness.
#' @param nmark Number of markers.
#' @param M Number of principal components.
#' @param ncovar Number of covariates.
#' @param lambda Additive predictor of time-varying survival covariates.
#' @param gamma Additive predictor of time-constant survival covariates.
#' @param alpha List of length nmark containing the additive predictors of the
#'   association.
#' @param mu List of length nmark containing the additive predictors of the 
#'   longitudinal part.
#' @param sigma Additive predictor of the variance.
#' @param tmax Maximal time point of observations.
#' @param seed Seed for reproducibility.
#' @param mfpc_args List containing the named arguments "type", "eFunType",
#'   "ignoreDeg", "eValType" of function simMultiFunData and "eValScale" for
#'   scaling the eigenvalues.
#' @param full Create a wide-format data.frame and 
simMultiJM <- function(nsub = 300, times = seq(0, 120, 1), probmiss = 0.75,
                       nmark = 2, M = 6, ncovar = 2,
                       lambda = function(t, x) {
                         1.4*log((t + 10)/1000) - 1.5
                       },
                       gamma = function(x) {
                         0.3*x[, 1]
                       },
                       alpha = rep(list(function(t, x) {
                         0.3 + 0*t
                       }), nmark),
                       mu = rep(list(function(t, x){
                         1.25 + 0.6*sin(x[, 2]) + (-0.01)*t
                       }), nmark),
                       sigma = function(t, x) {
                         0.3 + 0*t + I(x$marker == "m2")*0.2
                       }, 
                       tmax = NULL, seed = NULL, 
                       mfpc_args = list(type = "split", eFunType = "Poly",
                                        ignoreDeg = NULL, eValType = "linear",
                                        eValScale = 1),
                       full = FALSE, file = NULL){

  if(length(alpha) != length(mu)) {
    stop("alpha and mu must have same length.\n")
  }
  if(length(mu) != nmark) {
    stop("Predictors must be specified for all markers.\n")
  }
  if(is.null(tmax)){
    tmax <- max(times)
  }
  
  
  ## specify censoring function (aus CoxFlexBoost) changed into uniformly (as 
  ## too much censoring)
  ## added censoring at tmax
  censoring <- function(time, tmax){
    ## censoring times are independent uniformly distributed
    censor_time <- runif(n = length(time), min=0, max=1.5*tmax)
    censor_time <- ifelse(censor_time > tmax, tmax, censor_time)
    event <- (time <= censor_time)
    survtime <- apply(cbind(time, censor_time), 1, min)
    ## return matrix of observed survival times and event indicator
    return(cbind(survtime, event))
  }
  
  
  ## introduce random missings 
  ## (excluding the first measurement to ensure all subjects remain in the data)
  miss_fct <- function(data, prop, obstime = "obstime"){
    select <- which(data[[obstime]] > 0) 
    n <- length(select)
    n_miss <- round(prop*n, 0)
    miss <- sample(select, n_miss)
    data <- data[-miss,]
    return(data)
  }
  
  
  ## generate baseline covariates 
  gen_x <- function(nsub, ncovar){
    x <- matrix(data = NA, nrow = nsub, ncol = ncovar)
    for(i in seq_len(ncovar)) {
      x[, i] <- runif(n = nsub, min = -3, max = 3)
    }
    colnames(x) <- paste0("x", seq_len(ncovar))
    data.frame(x)
  } 
  
  
  ## Code from Clara Happ-Kurz' package funData and slightly adapted
  simMultiSplit <- function (argvals, M, eFunType, ignoreDeg = NULL, eValType, 
                             s) {
    if (any(c(length(M), length(eFunType), length(eValType)) != 1)) 
      stop("argvals, M, eFunType, eValType must all be of length 1!")
    p <- length(argvals)
    x <- vector("list", length = length(argvals))
    splitVals <- rep(NA, length(argvals) + 1)
    x[[1]] <- unlist(argvals[[1]])
    splitVals[1:2] <- c(0, length(x[[1]]))
    for (i in 2:p) {
      x[[i]] <- unlist(argvals[[i]])
      x[[i]] <- argvals[[i]] - min(argvals[[i]]) + max(x[[i - 1]])
      splitVals[i + 1] <- splitVals[i] + length(x[[i]])
    }
    f <- funData::eFun(unlist(x), M, ignoreDeg = ignoreDeg, type = eFunType)
    trueFuns <- vector("list", p)
    for (j in seq_len(p)) trueFuns[[j]] <- funData::funData(argvals[[j]], 
            s[j] * f@X[, (1 + splitVals[j]):splitVals[j + 1]])
    return(funData::multiFunData(trueFuns))
  }
  
  ## Code from Clara Happ-Kurz' package funData and slightly adapted
  simMultiWeight <- function (argvals, M, eFunType, ignoreDeg = NULL, eValType,
                              alpha) {
    p <- length(argvals)
    dimsSupp <- foreach::foreach(j = seq_len(p), .combine = "c") %do% {
      length(argvals[[j]])
    }
    if (any(dimsSupp > 2)) {
      stop(paste0("Function simMultiWeight: method is not implemented for ",
                  "objects of dimension > 2!"))
    }
    if (p > 1) {
      if (isTRUE(do.call(all.equal, lapply(M, prod)))) {
        Mtotal <- prod(M[[1]])
      }
      else stop("Function simMultiWeight: basis dimensions must be equal!")
    }
    else {
      Mtotal <- prod(M[[1]])
    }

    weight <- sqrt(alpha/sum(alpha))
    basis <- vector("list", p)
    for (j in seq_len(p)) {
      if (dimsSupp[j] == 1) 
        basis[[j]] <- weight[j] * funData::eFun(argvals[[j]][[1]], M = M[[j]], 
             ignoreDeg = ignoreDeg[[j]], type = eFunType[[j]])
      else basis[[j]] <- weight[j] * tensorProduct(
        funData::eFun(argvals[[j]][[1]], M = M[[j]][1], 
             ignoreDeg = ignoreDeg[[j]][[1]], type = eFunType[[j]][1]), 
              funData::eFun(argvals[[j]][[2]], M = M[[j]][2], 
                            ignoreDeg = ignoreDeg[[j]][[2]], 
                            type = eFunType[[j]][2]))
    }
    return(funData::multiFunData(basis))
  }
  
  
  ## generate the multivariate functional principal component basis
  mfpc <- function(argvals, mfpc_args, M) {
    switch(mfpc_args$type,
           split = simMultiSplit(argvals = argvals, M = M,
                                 eFunType = mfpc_args$eFunType,
                                 ignoreDeg = mfpc_args$ignoreDeg,
                                 eValType = mfpc_args$eValType,
                                 s = mfpc_args$mfpc_seed),
           weighted = simMultiWeight(argvals = argvals, M = M,
                                     eFunType = mfpc_args$eFunType,
                                     ignoreDeg = mfpc_args$ignoreDeg,
                                     eValType = mfpc_args$eValType,
                                     alpha = mfpc_args$mfpc_seed),
           stop(paste("Choose either 'split' or 'weighted' for the simulation",
                      "of multivariate functional data.")))
  }
  
  ## generate functional principal component based random effects
  gen_fpc <- function(times, nsub, M, mfpc_args, tmax, seed = NULL){
    if(!is.null(seed)) set.seed(seed)
    
    evals <- funData::eVal(M = M, type = mfpc_args$eValType)
    
    scores <- mvtnorm::rmvnorm(nsub, sigma = diag(mfpc_args$eValScale * evals),
                               method="chol")
    colnames(scores) <- paste0("s", 1:M)
    
    mfpc_seed <- switch(mfpc_args$type, 
                        "split" = sample(c(-1, 1), nmark, 0.5),
                        "weight" = stats::runif(nmark, 0.2, 0.8))
    b_set <- c(list(tmin = min(c(times, tmax)), tmax = tmax, M = M,
                    mfpc_seed = mfpc_seed, evals = evals), mfpc_args)
    return(list(scores, b_set))
  }
  
  ## compute predictors
  ## individual longitudinal trajectories
  mu_fun <-  function(time, x, r, mu, b_set){
    
    # duplicate scores for the multiple integration points
    if(is.null(dim(r))){
      r <- matrix(r, nrow = length(time), ncol = b_set$M, byrow=TRUE)
    }
    
    # Evaluate the functional principal component bases for different markers
    pc_bases <- lapply(mfpc(argvals = rep(list(c(b_set$tmin, time, b_set$tmax)),
                                          length(mu)),
                            mfpc_args = b_set, M = b_set$M),
                       function (fundat){
                         t(fundat@X)[-c(1, 2+length(time)), ]
                       })
    
    mapply(function (mu_k, pc_k) {
      mu_k(t = time, x = x) + apply(pc_k*r, 1, sum)
    }, mu_k = mu, pc_k = pc_bases, SIMPLIFY = FALSE)
  }
  
  
  
  ## full hazard (including nonlinear effect)
  hazard <-  function(time, x, r, ...){
    # mu_fun könnte auch hier als Objekt erstellt und unten ersetzt werden
    # ich brauche also nur eine Funktion, die mit gegeben time, x und r mu auf
    # allen markern auswertet
    exp(lambda(time, x) + gamma(x) + 
          Reduce('+', mapply(function(alpha_k, mu_k) {
            alpha_k(time, x)*mu_k
          }, alpha_k = alpha, mu_k = mu_fun(time, x, r, mu, b_set),
          SIMPLIFY = FALSE)))
  }
  
  
  # generate input
  id <- rep(1:nsub, each=length(times))
  if(!is.null(seed)){
    set.seed(seed)
  }
  x <- gen_x(nsub, ncovar = ncovar)
  

  temp <- gen_fpc(times = times, nsub = nsub, M = M, mfpc_args = mfpc_args,
                  tmax = tmax)
  r <- temp[[1]]
  b_set <- temp[[2]]
  
  data_short <- bamlss::rJM(hazard, censoring, x, r, tmin = times[1], tmax = tmax) 
  
  
  ## Create the full simulated data
  data_base <- cbind(id, data_short[id,], obstime = rep(times, nsub))
  i <- !duplicated(data_base$id)
  data_base$id <- as.factor(data_base$id)
  
  # gamma and lambda have only joint intercept which is estimated in
  # predictor gamma
  #????
  data_grid <- data.frame(survtime = times,
                          mu = seq(-0.5, 2.5, length.out = length(times)))
  f_lambda <- lambda(data_short$survtime)  
  f_gamma <- gamma(data_short[,grep("x[0-9]+", colnames(data_short))])
  data_base$lambda <- lambda(data_base$obstime) - mean(f_lambda)
  data_grid$lambda <- lambda(data_grid$survtime) - mean(f_lambda)
  data_base$gamma <- gamma(x[id, ]) + mean(f_lambda)
  data_grid$alpha <- cbind(data_grid,
                           alpha = do.call(cbind, 
                                           lapply(alpha, function(alpha_k) {
                                             alpha_k(data_grid$survtime, 0)
                                           })))
  #????
  data_long <- do.call(rbind, rep(list(data_base), nmark))
  data_long$marker <- factor(rep(paste0("m", seq_len(nmark)),
                                 each = length(id)))
  data_long$mu <- do.call(c, mu_fun(data_base$obstime, x[id, ], r[id,], mu,
                                    b_set))
  data_long$alpha <- do.call(c, lapply(alpha, function(alpha_k) {
                       alpha_k(data_base$obstime, x[id, ])
                     }))
  data_long$sigma <- sigma(t = data_long$obstime, x = data_long)
  fpcs <- do.call(rbind, lapply(mfpc(argvals = rep(list(data_base$obstime), 
                                                   nmark), mfpc_args = b_set,
                                     M = M), function (mark) {
                                       t(mark@X)
                                     }))
  data_long <- cbind(data_long,
                     fpc = fpcs,
                     wfpc = t(t(fpcs)*b_set$evals))
  
  # censoring                   
  data_long <- data_long[data_long$obstime <= data_long$survtime,]
  
  # saving data without longitudinal missings
  data_full <- data_long
  
  # inducing longitudinal missings
  data_long <- miss_fct(data_long, probmiss)
  
  # Draw longitudinal observations
  data_long$y <- rnorm(nrow(data_long), data_long$mu, sd = exp(data_long$sigma))
  
  
  if(full){
    d <- list(data=data_long, data_grid = data_grid, data_full = data_full)
  } else {
    d <- data_long
  }
  if(!is.null(file)) { 
    save(d, file = file) 
    invisible(d) 
  } else { 
    return(d) 
  } 
}
