

# New Simulation Function For Multivariate JMs Based On FPCs --------------


# Adapt the structure given by simJM function in bamlss

#' @param nsub Number of subjects.
#' @param times Vector of time points.
#' @param probmiss Probability of missingness.
#' @param ndim Number of dimensions.
#' @param M Number of principal components.
#' @param lambda Additive predictor of time-varying survival covariates.
#' @param gamma Additive predictor of time-constant survival covariates.
#' @param alpha List of length ndim containing the additive predictors of the
#'   association.
#' @param mu List of length ndim containing the additive predictors of the 
#'   longitudinal part.
#' @param sigma Additive predictor of the variance.
#' @param tmax Maximal time point of observations.
#' @param seed Seed for reproducibility.
#' @param mfpc_args List containing the named arguments "type", "eFunType",
#'   "ignoreDeg", and "eValType" of function simMultiFunData
simMultiJM <- function(nsub = 300, times = seq(0, 120, 1), probmiss = 0.75,
                       ndim = 2, M = 6,
                       lambda = function(t) {
                         1.4*log((time + 10)/1000) - 1.5
                       },
                       gamma = function(x1) {
                         0.3*x1
                       },
                       alpha = rep(list(function(t) {
                         3
                       }), ndim),
                       mu = rep(list(function(t, x){
                         1.25 + 0.6*sin(x) + (-0.01)*t
                       }), ndim),
                       sigma = function(t) {
                         0.3
                       }, 
                       tmax = NULL, seed = NULL, 
                       mfpc_args = list(type = "split", eFunType = "Poly",
                                        ignoreDeg = NULL, eValType = "linear"),
                  full = FALSE, file = NULL, nonlinear = FALSE, 
                  fac = FALSE, efun_type = "Poly", ignoreDeg = NULL, 
                  eval_type = "linear", eval_scale = 1){

  if(length(alpha) != length(mu)) {
    stop("alpha and mu must have same length.\n")
  }
  if(length(mu) != ndim) {
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
  ## (x1=survival covariate, x2=longitudinal covariate)
  gen_x <- function(nsub){
    x1 <- runif(nsub, -3, 3) 
    x2 <- runif(nsub, -3, 3)
    cbind(x1, x2)
  } 
  
  
  ## generate the multivariate functional principal component basis
  mfpc <- function(argvals, mfpc_args, M) {
    switch(mfpc_args$type,
           split = funData:::simMultiSplit(argvals = argvals, M = M,
                                           eFunType = mfpc_args$eFunType,
                                           ignoreDeg = mfpc_args$ignoreDeg,
                                           eValType = mfpc_args$eValType,
                                           N = 0),
           weighted = funData:::simMultiWeight(argvals = argvals, M = M,
                                               eFunType = mfpc_args$eFunType,
                                               ignoreDeg = mfpc_args$ignoreDeg,
                                               eValType = mfpc_args$eValType,
                                               N = 0),
           stop(paste("Choose either 'split' or 'weighted' for the simulation",
                      "of multivariate functional data."))
  }
  
  ## generate functional principal component based random effects
  gen_fpc <- function(times, nsub, M, eval_type = "Poly",
                      eval_scale = 1, tmax, efun_type = "linear", 
                      ignoreDeg = NULL, seed = NULL){
    if(!is.null(seed)) set.seed(seed)
    
    evals <- funData::eVal(M = M, type = eval_type)
    
    scores <- mvtnorm::rmvnorm(nsub, sigma = diag(eval_scale * evals),
                               method="chol")
    colnames(scores) <- paste0("s", 1:M)
    b_set <- list(tmin = min(c(times, tmax)),
                  tmax = tmax,
                  ignoreDeg = ignoreDeg,
                  type = efun_type)
    return(list(scores, b_set))
  }
  
  ## compute predictors
  ## individual longitudinal trajectories
  mu <-  function(time, x, r, M, b_set, long_setting){
    
    # allow to always extract coefficients as column
    if(is.null(dim(r))){
      r <- matrix(r, nrow=1)
    } 
    
    beta <- r[, -c(1:2)]
    # duplicate vector beta for the multiple integration points
    if(is.null(dim(beta))){
      beta <- matrix(beta, nrow = length(time), ncol = M, byrow=TRUE)
    }
    # TODO: Suppress warnings
    switch(long_setting,
     "linear" = (1.25 + r[, 1] + 0.6*sin(x) + (-0.01)*time + r[, 2]*0.02*time),
     "nonlinear" = (0.5 + r[, 1] + 0.6*sin(x) + 0.1*(time+2)*exp(-0.075*time)),
     "functional" = (0.5 + r[, 1] + 0.6*sin(x) + 0.1*(time+2)*exp(-0.075*time) + 
                       apply(splines::bs(time, M, b_set$knots,
                                         b_set$degree, b_set$intercept,
                                         b_set$Boundary.knots) * beta, 1, sum)),
     # Here adapt the longitudinal model
     "fpc" = (1.25 + 0.6*sin(x) + (-0.01)*time +
                apply(t(funData::eFun(argvals = c(b_set$tmin, time, b_set$tmax),
                   M = M, ignoreDeg = b_set$ignoreDeg,
                   type = b_set$type)@X)[-c(1, 2+length(time)), ]*beta, 1, sum)))
    
  }
  
  ## association between mu and log-hazard
  alpha <- function(time, alpha_setting){
    switch(alpha_setting,
           "zero" = 0*time,
           "constant" = 0*time + 1,
           "linear" = 1 - 0.015*time,
           "nonlinear" = cos((time-60)/20) + 1,
           "nonlinear2" = cos((time-33)/33))
  }
  
  alpha_nonlin <- function(time, alpha_setting, x2, x3, r, M, b_set, long_setting){
    switch(alpha_setting,
           "zero" = 0*time,
           "constant" = 0*time + 1,
           "linear" = 1 * mu(time, x2, r, M, b_set, long_setting),
           "nonlinear" = -0.1 * (mu(time, x2, r, M, b_set, long_setting) + 3)^2 +
             mu(time, x2, r, M, b_set, long_setting) + 1.8,
           "nonlinear2" = x3 * (-0.1 * (mu(time, x2, r, M, b_set, long_setting) + 3)^2 +
                                  mu(time, x2, r, M, b_set, long_setting) + 1.8) + 
             (1 - x3) * (0.1 * (mu(time, x2, r, M, b_set, long_setting) - 3)^2 +
                           0.75 * mu(time, x2, r, M, b_set, long_setting) - 0.8))
  }
  
  # used only for setting up final data
  alpha_nonlin_simple <- function(alpha_setting, x3, mu){
    switch(alpha_setting,
           "zero" = 0,
           "constant" = 1,
           "linear" = 1 * mu,
           "nonlinear" = -0.1 * (mu + 3)^2 + mu + 1.8,
           "nonlinear2" = x3 * (-0.1 * (mu + 3)^2 +  mu + 1.8) + 
             (1 - x3) * (0.1 * (mu - 3)^2 + 0.75 * mu - 0.8))
  }
  
  
  ## baseline hazard
  lambda <- pred_lambda[[1]]
  
  ## baseline covariate
  gamma <-  pred_gamma[[1]]
  
  ## full hazard (including nonlinear effect)
  hazard <-  function(time, x, r, ...){
    if(nonlinear){
      exp(lambda(time) + gamma(x[1], nonlinear) + 
            alpha_nonlin(time, alpha_setting, x[2], x[3], r, M, b_set, long_setting))
    } else{
      exp(lambda(time) + gamma(x[1]) + 
            alpha(time, alpha_setting)*mu(time, x[2], r, M, b_set, long_setting))
    }
  }
  
  
  # generate input
  id <- rep(1:nsub, each=length(times))
  if(!is.null(seed)){
    set.seed(seed)
  }
  r <- gen_r(nsub)
  x <- gen_x(nsub)
  if(!fac) x[, 3] <- rep(1, nsub)
  
  if (long_setting == "fpc") {
    temp <- gen_fpc(times = times, nsub = nsub, M = M, 
                    eval_type = eval_type, eval_scale = eval_scale,
                    tmax = tmax, efun_type = efun_type, ignoreDeg = ignoreDeg)
  } else {
    temp <- gen_b(times, nsub, M=M, l=c(1,5))
  }
  r <- cbind(r, temp[[1]])
  b_set <- temp[[2]]
  
  data_short <- rJM(hazard, censoring, x, r, tmin = times[1], tmax = tmax) 
  
  data_long <- cbind(id, data_short[id,], obstime = rep(times, nsub))
  data_grid <- data.frame(survtime = times,
                          mu = seq(-0.5, 2.5, length.out = length(times)))
  
  i <- !duplicated(data_long$id)
  
  
  # gamma and lambda have only joint intercept which is estimated in predictor gamma
  data_long$mu <- mu(data_long$obstime, data_long$x2, r[id,], M, b_set, long_setting)
  f_lambda <- lambda(data_long$survtime)[i]   
  f_gamma <- gamma(data_long$x1, nonlinear)[i]
  data_long$lambda <- lambda(data_long$survtime) - mean(f_lambda)
  data_grid$lambda <- lambda(data_grid$survtime) - mean(f_lambda)
  if(nonlinear){
    data_long$gamma <- gamma(data_long$x1, nonlinear) + mean(f_lambda)
    data_long$surv_mu <- mu(data_long$survtime, data_long$x2, r[id,], M, b_set, long_setting)
    data_long$alpha <- alpha_nonlin_simple(alpha_setting, data_long$x3, data_long$surv_mu)
    data_long$alpha_l <- alpha_nonlin_simple(alpha_setting, data_long$x3, data_long$mu)
    data_grid$mu <- seq(-0.5, 2.5, length.out = nrow(data_grid))
    if(fac){
      data_grid$alpha1 <- alpha_nonlin_simple(alpha_setting, rep(1, nrow(data_grid)), data_grid$mu)
      data_grid$alpha0 <- alpha_nonlin_simple(alpha_setting, rep(0, nrow(data_grid)), data_grid$mu)
    } else {
      data_grid$alpha <- alpha_nonlin_simple(alpha_setting, rep(1, nrow(data_grid)), data_grid$mu)
    }
    
  } else {
    data_long$alpha <- alpha(data_long$survtime, alpha_setting)
    data_grid$alpha <- alpha(data_grid$survtime, alpha_setting)
    data_long$gamma <- gamma(data_long$x1, nonlinear) + mean(f_lambda)
  }
  
  data_long$dmu <- dmu(data_long$obstime, r[id,], M, b_set, long_setting)
  data_long$id <- as.factor(data_long$id)
  data_long$sigma <- rep(log(sigma), nrow(data_long))
  data_long$x3 <- as.factor(data_long$x3)
  
  # censoring                   
  data_long <- data_long[data_long$obstime <= data_long$survtime,]
  
  # saving data without longitudinal missings
  data_full <- data_long
  
  # inducing longitudinal missings
  data_long <- miss_fct(data_long, probmiss)
  
  # Draw longitudinal observations
  data_long$y <- rnorm(nrow(data_long), data_long$mu, sigma)
  
  ygrid <- quantile(data_long$y, probs = seq(0.025, 0.975, 0.025))
  # adjust predictions with constraint median(y)
  if(nonlinear){
    # constraint acts on group 0. to be checked
    alpha_constraint <- mean(alpha_nonlin_simple(alpha_setting, 0, ygrid))
    data_long$alpha <- data_long$alpha - alpha_constraint
    data_long$alpha_l <- data_long$alpha_l - alpha_constraint
    if(fac){
      data_grid$alpha1 <- data_grid$alpha1 - alpha_constraint
      data_grid$alpha0 <- data_grid$alpha0 - alpha_constraint
    } else {
      data_grid$alpha <- data_grid$alpha - alpha_constraint
    }
    data_long$gamma <- data_long$gamma + alpha_constraint
    data_long$alpha_constraint <- alpha_constraint
  } 
  
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
