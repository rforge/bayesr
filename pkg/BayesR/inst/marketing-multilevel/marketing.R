library("foreign")
library("BayesR")


dir <- path.expand("~/diss/marketing")
datafiles <- file.path(dir, grep(".dta", list.files(dir), value = TRUE))


## hypident = TRUE, contr = 2
for(file in datafiles) {
  dat <- read.dta(file)
  bsn <- paste("brand-", gsub("_cv.dta", "", basename(file)), "-hypident-constrained", sep = "")

  b <- gibbs(labsatz ~ m(outlet, s(preis_own, bs = "ps", k = 20, xt = list(constr = 2)), hypident = TRUE), 
    iter = 12000, burnin = 2000, thinning = 10, data = dat)

  dic <- b$DIC
  k <- 1L
  fits <- NULL
  for(id in unique(dat$outlet)) {
    check <- dat$outlet == id
    response <- dat$labsatz[check][order(dat$preis_own[check])]
    fits <- rbind(fits, cbind(response, id, b$fout[[1L]][[k]][, 1L:2L]))
    k <- k + 1L
  }
  colnames(fits) <- c("labsatz", "outlet", "preis_own", "fitted")
  fits <- as.data.frame(fits)
  fits$fitted <- fits$fitted + coef(b)[1L, 1L]

  save(b, fits, dic, file = file.path(dir, paste(bsn, ".rda", sep = "")), compress = "xz")
  write.table(fits, file = file.path(dir, paste(bsn, "-fitted.raw", sep = "")),
    row.names = FALSE, quote = FALSE)
  write.table(dic, file = file.path(dir, paste(bsn, "-dic.raw", sep = "")),
    row.names = FALSE, quote = FALSE)
}


## hypident = FALSE, contr = 2
for(file in datafiles) {
  dat <- read.dta(file)
  bsn <- paste("brand-", gsub("_cv.dta", "", basename(file)), "-nohypident-constrained", sep = "")

  b <- gibbs(labsatz ~ m(outlet, s(preis_own, bs = "ps", k = 20, xt = list(constr = 2)), hypident = FALSE), 
    iter = 12000, burnin = 2000, thinning = 10, data = dat)

  dic <- b$DIC
  k <- 1L
  fits <- NULL
  for(id in unique(dat$outlet)) {
    check <- dat$outlet == id
    response <- dat$labsatz[check][order(dat$preis_own[check])]
    fits <- rbind(fits, cbind(response, id, b$fout[[1L]][[k]][, 1L:2L]))
    k <- k + 1L
  }
  colnames(fits) <- c("labsatz", "outlet", "preis_own", "fitted")
  fits <- as.data.frame(fits)
  fits$fitted <- fits$fitted + coef(b)[1L, 1L]

  save(b, fits, dic, file = file.path(dir, paste(bsn, ".rda", sep = "")), compress = "xz")
  write.table(fits, file = file.path(dir, paste(bsn, "-fitted.raw", sep = "")),
    row.names = FALSE, quote = FALSE)
  write.table(dic, file = file.path(dir, paste(bsn, "-dic.raw", sep = "")),
    row.names = FALSE, quote = FALSE)
}


## hypident = TRUE, unconstrained
for(file in datafiles) {
  dat <- read.dta(file)
  bsn <- paste("brand-", gsub("_cv.dta", "", basename(file)), "-hypident-unconstrained", sep = "")

  b <- gibbs(labsatz ~ m(outlet, s(preis_own, bs = "ps", k = 20), hypident = TRUE), 
    iter = 12000, burnin = 2000, thinning = 10, data = dat)

  dic <- b$DIC
  k <- 1L
  fits <- NULL
  for(id in unique(dat$outlet)) {
    check <- dat$outlet == id
    response <- dat$labsatz[check][order(dat$preis_own[check])]
    fits <- rbind(fits, cbind(response, id, b$fout[[1L]][[k]][, 1L:2L]))
    k <- k + 1L
  }
  colnames(fits) <- c("labsatz", "outlet", "preis_own", "fitted")
  fits <- as.data.frame(fits)
  fits$fitted <- fits$fitted + coef(b)[1L, 1L]

  save(b, fits, dic, file = file.path(dir, paste(bsn, ".rda", sep = "")), compress = "xz")
  write.table(fits, file = file.path(dir, paste(bsn, "-fitted.raw", sep = "")),
    row.names = FALSE, quote = FALSE)
  write.table(dic, file = file.path(dir, paste(bsn, "-dic.raw", sep = "")),
    row.names = FALSE, quote = FALSE)
}


## hypident = FALSE, unconstrained
for(file in datafiles) {
  dat <- read.dta(file)
  bsn <- paste("brand-", gsub("_cv.dta", "", basename(file)), "-nohypident-unconstrained", sep = "")

  b <- gibbs(labsatz ~ m(outlet, s(preis_own, bs = "ps", k = 20), hypident = FALSE), 
    iter = 12000, burnin = 2000, thinning = 10, data = dat)

  dic <- b$DIC
  k <- 1L
  fits <- NULL
  for(id in unique(dat$outlet)) {
    check <- dat$outlet == id
    response <- dat$labsatz[check][order(dat$preis_own[check])]
    fits <- rbind(fits, cbind(response, id, b$fout[[1L]][[k]][, 1L:2L]))
    k <- k + 1L
  }
  colnames(fits) <- c("labsatz", "outlet", "preis_own", "fitted")
  fits <- as.data.frame(fits)
  fits$fitted <- fits$fitted + coef(b)[1L, 1L]

  save(b, fits, dic, file = file.path(dir, paste(bsn, ".rda", sep = "")), compress = "xz")
  write.table(fits, file = file.path(dir, paste(bsn, "-fitted.raw", sep = "")),
    row.names = FALSE, quote = FALSE)
  write.table(dic, file = file.path(dir, paste(bsn, "-dic.raw", sep = "")),
    row.names = FALSE, quote = FALSE)
}


#rm(list = ls())
#dir <- path.expand("~/diss/marketing")
#fitfiles <- file.path(dir, grep(".rda", list.files(dir), value = TRUE))
#load(fitfiles[8L])

#for(id in unique(fits$outlet)) {
#  plot(fits$preis_own, fits$labsatz, 
#    ylim = range(c(fits$labsatz, fits$fitted)), 
#    type = "n", xlab = "preis_own", ylab = "labsatz")
#  check <- id == fits$outlet
#  points(fits$preis_own[check], fits$labsatz[check])
#  lines(fits$fitted[check] ~ fits$preis_own[check], col = "green3")
#  Sys.sleep(0.7)
#}
