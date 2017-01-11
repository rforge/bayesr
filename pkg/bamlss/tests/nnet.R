library("bamlss")


derivs_nn <- function(k = 1)
{
  fun <- "w1z * x"
  j <- 2
  for(i in 1:k) {
    fun <- c(fun, paste("w", j, "z * 1 / (1 + exp(-(w", j + 1, "z + w", j + 2, "z * x)))", sep = ""))
    j <- j + 3
  }
  w <- paste("w", 1:(j - 1), "z", sep = "")
  w2 <- paste("w", 1:(j - 1), sep = "")
  fun <- paste(fun, collapse = " + ")
  efun <- parse(text = fun)
  grad <- hess <- list()
  for(j in seq_along(w)) {
    grad[[w2[j]]] <- D(efun, w[j])
    for(i in seq_along(w)) {
      if(i <= j)
        hess[[paste(w2[j], ",", w2[i], sep = "")]] <- D(grad[[w2[j]]], w[i])
    }
  }
  for(j in seq_along(grad)) {
    grad[[j]] <- paste(deparse(grad[[j]]), collapse = "")
    for(i in seq_along(w))
      grad[[j]] <- gsub(paste("w", i, "z", sep = ""), paste("w[", i, "]", sep = ""), grad[[j]])
    grad[[j]] <- eval(parse(text = paste("function(w) {", grad[[j]], "}")))
  }
  for(j in seq_along(hess)) {
    hess[[j]] <- paste(deparse(hess[[j]]), collapse = "")
    for(i in seq_along(w)) {
      hess[[j]] <- gsub(paste("w", i, "z", sep = ""), paste("w[", i, "]", sep = ""), hess[[j]])
    }
    hess[[j]] <- eval(parse(text = paste("function(w) {", hess[[j]], "}")))
  }
  list("grad" = grad, "hess" = hess)
}


fit.nn <- function(x, y, k = 10, start = NULL,
  method = c("optim", "nr"), plot = FALSE,
  lambda = .Machine$double.eps^0.5, control = list(...), ...)
{
  method <- match.arg(method)

  if(is.null(control$maxit))
    control$maxit <- 2000
  if(is.null(control$reltol))
    control$reltol <- 0.00001 #sqrt(.Machine$double.eps)
  if(is.null(control$verbose))
    control$verbose <- TRUE

  if(!is.null(start))
    k <- length(start)

  gnn <- derivs_nn(k = k)$grad

  fit_nn <- function(w) {
    w[is.na(w)] <- .Machine$double.eps
    eta <- w[1] * x
    j <- 2
    for(i in 1:k) {
      eta <- eta + w[j] * 1 / (1 + exp(-(w[j + 1] + w[j + 2] * x)))
      j <- j + 3
    }
    eta
  }

  objfun_nn <- function(w) {
    w[is.na(w)] <- .Machine$double.eps
    fit <- fit_nn(w)
    if(plot) {
      plot(x, y, ylim = range(c(y, fit)))
      ox <- order(x)
      lines(fit[ox] ~ x[ox], lwd = 2)
    }
    sum((y - fit)^2) + lambda * sum(w^2)
  }

  grad_nn <- function(w) {
    w[is.na(w)] <- .Machine$double.eps
    dfdeta <- -(2 * (y - fit_nn(w)))
    grad <- rep(0, length(w))
    for(i in seq_along(w))
      grad[i] <- sum(gnn[[i]](w) * dfdeta)
    grad + lambda * w
  }

  hess_nn <- function(w) {
    w[is.na(w)] <- .Machine$double.eps
    m <- length(w)
    h <- matrix(0, m, m)
    for(i in 1:m) {
      for(j in 1:m) {
        if(j <= i) {
          h[i, j] <- crossprod(gnn[[i]](w), gnn[[j]](w))
        }
        h[j, i] <- h[i, j]
      }
    }
    h + lambda * diag(length(w))
  }

  w <- if(!is.null(start)) start else runif(k * 3 + 1, 0, 1)

  rval <- list()
  if(method == "optim") {
    control$verbose <- NULL
    opt <- optim(w, fn = objfun_nn, gr = grad_nn, method = "BFGS", control = control, hessian = TRUE)
    rval$coefficients <- opt$par
    rval$hessian <- opt$hessian
    rval$fitted.values <- fit_nn(opt$par)
    rval$converged <- opt$convergence < 1
  } else {
    eps <- control$reltol + 1
    iter <- 0
    fit <- fit_nn(w)
    while(eps > control$reltol & iter < control$maxit) {
      fit0 <- fit
      g <- grad_nn(w)
      h <- chol2inv(chol(hess_nn(w) + diag(.Machine$double.eps^0.5, length(w))))
      foo <- function(nu) {
        sum((y - fit_nn(drop(w - nu * h %*% g)))^2)
      }
      nu <- optimize(foo, lower = 0, upper = 1)$minimum
      w <- w - drop(nu * h %*% g)
      fit <- fit_nn(w)
      eps <- mean(abs((fit - fit0) / fit0))
      iter <- iter + 1
      if(control$verbose) {
        cat("\r")
        vtxt <- paste(
          " eps ", bamlss:::fmt(eps, width = 6, digits = 6),
          " iteration ", formatC(iter, width = nchar(control$maxit)), sep = "")
        cat(vtxt)
      }
      if(plot) {
        plot(x, y, ylim = range(c(y, fit)))
        ox <- order(x)
        lines(fit[ox] ~ x[ox], lwd = 2)
      }
    }
    if(control$verbose) cat("\n")
    rval$coefficients <- w
    rval$hessian <- hess_nn(w) + diag(.Machine$double.eps^0.5, length(w))
    rval$fitted.values <- fit_nn(w)
  }

  rval
}


f <- bamlss:::simfun("complicated")

n <- 300
x <- seq(0, 1, length = n)
y <- 1.2 + f(x) + rnorm(n, sd = 0.1)

b <- fit.nn(x, y, k = 10, method = "optim", plot = TRUE)

