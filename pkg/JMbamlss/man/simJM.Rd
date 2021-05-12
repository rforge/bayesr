\name{simJM}
\alias{simJM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simJM(nsub = 300, times = seq(0, 120, 1), probmiss = 0.75, long_setting = "functional", alpha_setting = if (nonlinear) "linear" else "nonlinear", dalpha_setting = "zero", sigma = 0.3, long_df = 6, tmax = NULL, seed = NULL, full = FALSE, file = NULL, nonlinear = FALSE, fac = FALSE, efun_type = "Poly", ignoreDeg = NULL, eval_type = "linear", eval_scale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsub}{
%%     ~~Describe \code{nsub} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{probmiss}{
%%     ~~Describe \code{probmiss} here~~
}
  \item{long_setting}{
%%     ~~Describe \code{long_setting} here~~
}
  \item{alpha_setting}{
%%     ~~Describe \code{alpha_setting} here~~
}
  \item{dalpha_setting}{
%%     ~~Describe \code{dalpha_setting} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{long_df}{
%%     ~~Describe \code{long_df} here~~
}
  \item{tmax}{
%%     ~~Describe \code{tmax} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{full}{
%%     ~~Describe \code{full} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{nonlinear}{
%%     ~~Describe \code{nonlinear} here~~
}
  \item{fac}{
%%     ~~Describe \code{fac} here~~
}
  \item{efun_type}{
%%     ~~Describe \code{efun_type} here~~
}
  \item{ignoreDeg}{
%%     ~~Describe \code{ignoreDeg} here~~
}
  \item{eval_type}{
%%     ~~Describe \code{eval_type} here~~
}
  \item{eval_scale}{
%%     ~~Describe \code{eval_scale} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (nsub = 300, times = seq(0, 120, 1), probmiss = 0.75, 
    long_setting = "functional", alpha_setting = if (nonlinear) "linear" else "nonlinear", 
    dalpha_setting = "zero", sigma = 0.3, long_df = 6, tmax = NULL, 
    seed = NULL, full = FALSE, file = NULL, nonlinear = FALSE, 
    fac = FALSE, efun_type = "Poly", ignoreDeg = NULL, eval_type = "linear", 
    eval_scale = 1) 
{
    if (is.null(tmax)) {
        tmax <- max(times)
    }
    if (nonlinear & alpha_setting == "nonlinear2") 
        fac <- TRUE
    if (nonlinear & fac & alpha_setting == "nonlinear") 
        alpha_setting <- "nonlinear2"
    censoring <- function(time, tmax) {
        censor_time <- runif(n = length(time), min = 0, max = 1.5 * 
            tmax)
        censor_time <- ifelse(censor_time > tmax, tmax, censor_time)
        event <- (time <= censor_time)
        survtime <- apply(cbind(time, censor_time), 1, min)
        return(cbind(survtime, event))
    }
    miss_fct <- function(data, prop, obstime = "obstime") {
        select <- which(data[[obstime]] > 0)
        n <- length(select)
        n_miss <- round(prop * n, 0)
        miss <- sample(select, n_miss)
        data <- data[-miss, ]
        return(data)
    }
    gen_x <- function(nsub) {
        x1 <- runif(nsub, -3, 3)
        x2 <- runif(nsub, -3, 3)
        x3 <- rbinom(nsub, 1, 0.5)
        cbind(x1, x2, x3)
    }
    gen_r <- function(nsub) {
        r1 <- r2 <- matrix(NA, nrow = nsub, ncol = 2)
        r1 <- rnorm(nsub, 0, 0.25)
        r2 <- rnorm(nsub, 0, 0.4)
        cbind(r1, r2)
    }
    gen_b <- function(times, nsub, long_df, pen = 2, l = c(1, 
        1), seed = NULL) {
        if (!is.null(seed)) 
            set.seed(seed)
        makeDiffOp <- function(degree, dim) {
            if (degree == 0) {
                return(diag(dim))
            }
            else {
                return(diff(makeDiffOp(degree - 1, dim)))
            }
        }
        Pi <- l[1] * diag(long_df)
        Pt <- l[2] * crossprod(makeDiffOp(pen[1], long_df))
        P <- 0.1 * diag(long_df) + Pt + Pi
        coef <- mvtnorm::rmvnorm(nsub, sigma = solve(P), method = "chol")
        colnames(coef) <- paste0("b", 1:long_df)
        bt <- splines::bs(times, df = long_df, intercept = FALSE)
        b_set <- list(knots = attr(bt, "knots"), Boundary.knots = attr(bt, 
            "Boundary.knots"), degree = attr(bt, "degree"), intercept = attr(bt, 
            "intercept"))
        return(list(coef, b_set))
    }
    dbs <- function(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE, 
        Boundary.knots = range(x), eps = 1e-07) {
        ex <- pmax(abs(x), 1)
        x1 <- x + eps * ex
        x2 <- x - eps * ex
        bs.xeps1 <- suppressWarnings(splines::bs(x1, df, knots, 
            degree, intercept, Boundary.knots))
        bs.xeps2 <- suppressWarnings(splines::bs(x2, df, knots, 
            degree, intercept, Boundary.knots))
        out <- (bs.xeps1 - bs.xeps2)/c(x1 - x2)
        out
    }
    gen_fpc <- function(times, nsub, long_df, eval_type = "Poly", 
        eval_scale = 1, tmax, efun_type = "linear", ignoreDeg = NULL, 
        seed = NULL) {
        if (!is.null(seed)) 
            set.seed(seed)
        evals <- funData::eVal(M = long_df, type = eval_type)
        scores <- mvtnorm::rmvnorm(nsub, sigma = diag(eval_scale * 
            evals), method = "chol")
        colnames(scores) <- paste0("s", 1:long_df)
        b_set <- list(tmin = min(c(times, tmax)), tmax = tmax, 
            ignoreDeg = ignoreDeg, type = efun_type)
        return(list(scores, b_set))
    }
    mu <- function(time, x, r, long_df, b_set, long_setting) {
        if (is.null(dim(r))) {
            r <- matrix(r, nrow = 1)
        }
        beta <- r[, -c(1:2)]
        if (is.null(dim(beta))) {
            beta <- matrix(beta, nrow = length(time), ncol = long_df, 
                byrow = TRUE)
        }
        switch(long_setting, linear = (1.25 + r[, 1] + 0.6 * 
            sin(x) + (-0.01) * time + r[, 2] * 0.02 * time), 
            nonlinear = (0.5 + r[, 1] + 0.6 * sin(x) + 0.1 * 
                (time + 2) * exp(-0.075 * time)), functional = (0.5 + 
                r[, 1] + 0.6 * sin(x) + 0.1 * (time + 2) * exp(-0.075 * 
                time) + apply(splines::bs(time, long_df, b_set$knots, 
                b_set$degree, b_set$intercept, b_set$Boundary.knots) * 
                beta, 1, sum)), fpc = (1.25 + 0.6 * sin(x) + 
                (-0.01) * time + apply(t(funData::eFun(argvals = c(b_set$tmin, 
                time, b_set$tmax), M = long_df, ignoreDeg = b_set$ignoreDeg, 
                type = b_set$type)@X)[-c(1, 2 + length(time)), 
                ] * beta, 1, sum)))
    }
    dmu <- function(time, r, long_df, b_set, long_setting) {
        if (is.null(dim(r))) {
            r <- matrix(r, nrow = 1)
        }
        beta <- r[, -c(1:2)]
        if (is.null(dim(beta))) {
            beta <- matrix(beta, nrow = length(time), ncol = long_df, 
                byrow = TRUE)
        }
        switch(long_setting, fpc = , simple = (-0.02) + 0 * time, 
            linear = (-0.01 + r[, 2] * 0.02) + 0 * time, nonlinear = (0.085 - 
                0.0075 * time) * exp(-0.075 * time), functional = (0.085 - 
                0.0075 * time) * exp(-0.075 * time) + apply(dbs(time, 
                long_df, b_set$knots, b_set$degree, b_set$intercept, 
                b_set$Boundary.knots) * beta, 1, sum))
    }
    alpha <- function(time, alpha_setting) {
        switch(alpha_setting, zero = 0 * time, constant = 0 * 
            time + 1, linear = 1 - 0.015 * time, nonlinear = cos((time - 
            60)/20) + 1, nonlinear2 = cos((time - 33)/33))
    }
    alpha_nonlin <- function(time, alpha_setting, x2, x3, r, 
        long_df, b_set, long_setting) {
        switch(alpha_setting, zero = 0 * time, constant = 0 * 
            time + 1, linear = 1 * mu(time, x2, r, long_df, b_set, 
            long_setting), nonlinear = -0.1 * (mu(time, x2, r, 
            long_df, b_set, long_setting) + 3)^2 + mu(time, x2, 
            r, long_df, b_set, long_setting) + 1.8, nonlinear2 = x3 * 
            (-0.1 * (mu(time, x2, r, long_df, b_set, long_setting) + 
                3)^2 + mu(time, x2, r, long_df, b_set, long_setting) + 
                1.8) + (1 - x3) * (0.1 * (mu(time, x2, r, long_df, 
            b_set, long_setting) - 3)^2 + 0.75 * mu(time, x2, 
            r, long_df, b_set, long_setting) - 0.8))
    }
    alpha_nonlin_simple <- function(alpha_setting, x3, mu) {
        switch(alpha_setting, zero = 0, constant = 1, linear = 1 * 
            mu, nonlinear = -0.1 * (mu + 3)^2 + mu + 1.8, nonlinear2 = x3 * 
            (-0.1 * (mu + 3)^2 + mu + 1.8) + (1 - x3) * (0.1 * 
            (mu - 3)^2 + 0.75 * mu - 0.8))
    }
    dalpha <- function(time, dalpha_setting) {
        switch(dalpha_setting, zero = 0 * time, constant = 0 * 
            time + 10, linear = 6 - 0.015 * time, nonlinear = 50 + 
            55 * sin((time)/20), nonlinear2 = 50 + 55 * sin((time)/20))
    }
    lambda <- function(time) {
        1.4 * log((time + 10)/1000) - 1.5
    }
    gamma <- function(x1, nonlinear = TRUE) {
        if (nonlinear) {
            0.3 * x1
        }
        else {
            sin(x1)
        }
    }
    hazard <- function(time, x, r, ...) {
        if (nonlinear) {
            exp(lambda(time) + gamma(x[1], nonlinear) + alpha_nonlin(time, 
                alpha_setting, x[2], x[3], r, long_df, b_set, 
                long_setting) + dalpha(time, dalpha_setting) * 
                dmu(time, r, long_df, b_set, long_setting))
        }
        else {
            exp(lambda(time) + gamma(x[1]) + alpha(time, alpha_setting) * 
                mu(time, x[2], r, long_df, b_set, long_setting) + 
                dalpha(time, dalpha_setting) * dmu(time, r, long_df, 
                  b_set, long_setting))
        }
    }
    id <- rep(1:nsub, each = length(times))
    if (!is.null(seed)) {
        set.seed(seed)
    }
    r <- gen_r(nsub)
    x <- gen_x(nsub)
    if (!fac) 
        x[, 3] <- rep(1, nsub)
    if (long_setting == "fpc") {
        temp <- gen_fpc(times = times, nsub = nsub, long_df = long_df, 
            eval_type = eval_type, eval_scale = eval_scale, tmax = tmax, 
            efun_type = efun_type, ignoreDeg = ignoreDeg)
    }
    else {
        temp <- gen_b(times, nsub, long_df = long_df, l = c(1, 
            5))
    }
    r <- cbind(r, temp[[1]])
    b_set <- temp[[2]]
    data_short <- rJM(hazard, censoring, x, r, tmin = times[1], 
        tmax = tmax)
    data_long <- cbind(id, data_short[id, ], obstime = rep(times, 
        nsub))
    data_grid <- data.frame(survtime = times, mu = seq(-0.5, 
        2.5, length.out = length(times)))
    i <- !duplicated(data_long$id)
    data_long$dalpha <- dalpha(data_long$survtime, dalpha_setting)
    data_grid$dalpha <- dalpha(data_grid$survtime, dalpha_setting)
    data_long$mu <- mu(data_long$obstime, data_long$x2, r[id, 
        ], long_df, b_set, long_setting)
    f_lambda <- lambda(data_long$survtime)[i]
    f_gamma <- gamma(data_long$x1, nonlinear)[i]
    data_long$lambda <- lambda(data_long$survtime) - mean(f_lambda)
    data_grid$lambda <- lambda(data_grid$survtime) - mean(f_lambda)
    if (nonlinear) {
        data_long$gamma <- gamma(data_long$x1, nonlinear) + mean(f_lambda)
        data_long$surv_mu <- mu(data_long$survtime, data_long$x2, 
            r[id, ], long_df, b_set, long_setting)
        data_long$alpha <- alpha_nonlin_simple(alpha_setting, 
            data_long$x3, data_long$surv_mu)
        data_long$alpha_l <- alpha_nonlin_simple(alpha_setting, 
            data_long$x3, data_long$mu)
        data_grid$mu <- seq(-0.5, 2.5, length.out = nrow(data_grid))
        if (fac) {
            data_grid$alpha1 <- alpha_nonlin_simple(alpha_setting, 
                rep(1, nrow(data_grid)), data_grid$mu)
            data_grid$alpha0 <- alpha_nonlin_simple(alpha_setting, 
                rep(0, nrow(data_grid)), data_grid$mu)
        }
        else {
            data_grid$alpha <- alpha_nonlin_simple(alpha_setting, 
                rep(1, nrow(data_grid)), data_grid$mu)
        }
    }
    else {
        data_long$alpha <- alpha(data_long$survtime, alpha_setting)
        data_grid$alpha <- alpha(data_grid$survtime, alpha_setting)
        data_long$gamma <- gamma(data_long$x1, nonlinear) + mean(f_lambda)
    }
    data_long$dmu <- dmu(data_long$obstime, r[id, ], long_df, 
        b_set, long_setting)
    data_long$id <- as.factor(data_long$id)
    data_long$sigma <- rep(log(sigma), nrow(data_long))
    data_long$x3 <- as.factor(data_long$x3)
    data_long <- data_long[data_long$obstime <= data_long$survtime, 
        ]
    data_full <- data_long
    data_long <- miss_fct(data_long, probmiss)
    data_long$y <- rnorm(nrow(data_long), data_long$mu, sigma)
    ygrid <- quantile(data_long$y, probs = seq(0.025, 0.975, 
        0.025))
    if (nonlinear) {
        alpha_constraint <- mean(alpha_nonlin_simple(alpha_setting, 
            0, ygrid))
        data_long$alpha <- data_long$alpha - alpha_constraint
        data_long$alpha_l <- data_long$alpha_l - alpha_constraint
        if (fac) {
            data_grid$alpha1 <- data_grid$alpha1 - alpha_constraint
            data_grid$alpha0 <- data_grid$alpha0 - alpha_constraint
        }
        else {
            data_grid$alpha <- data_grid$alpha - alpha_constraint
        }
        data_long$gamma <- data_long$gamma + alpha_constraint
        data_long$alpha_constraint <- alpha_constraint
    }
    if (full) {
        d <- list(data = data_long, data_grid = data_grid, data_full = data_full)
    }
    else {
        d <- data_long
    }
    if (!is.null(file)) {
        save(d, file = file)
        invisible(d)
    }
    else {
        return(d)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
