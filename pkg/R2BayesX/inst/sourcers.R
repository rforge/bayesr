dir <- "/home/c403129/svn/bayesr/pkg/R2BayesX/R"
## dir <- "J:/c403/stat/R2BayesX/R"
invisible(sapply(paste(dir, "/", list.files(dir), sep = ""), source))


names(columb.polys) <- 1:length(columb.polys)
plotmap(map = columb.polys, x = columb$crime, id = columb$district, values = T)
cbind(columb$crime, f2int(columb$district))

library("R2BayesX")
## Load Columbus Ohio crime data 
## and polygon shape file
data("columb")     
data("columb.polys")

## need to adapt pylgon names
## for comparison
columb$district <- with(columb, x2int(district))
names(columb.polys) <- x2int(names(columb.polys))

## neighbourhood structure info for MRF
xt <- list(polys = columb.polys) 

## estimate models with
## mgcv REML and BayesX MCMC 
b1 <- gam(crime ~ s(district, bs = "mrf", xt = xt), data = columb, method = "REML")
b2 <- bayesx(crime ~ s(district, bs = "mrf", xt = xt), data = columb, method = "REML")

fit.b2 <- fitted(b2, term = "s(district)")[[1L]][,1L:2L] 
fit.b2[,2L] <- fit.b2[,2L] + coef(b2)[1L, 1L]

plotmap(map = columb.polys, x = predict(b1, terms = 1), symmetric = FALSE, values = TRUE)
plotmap(map = columb.polys, x = fit.b2, symmetric = FALSE, values = TRUE)


set.seed(333)
     
## simulate some geographical data
data("MunichBnd")
N <- length(MunichBnd); n <- N*5
     
## regressors
dat <- data.frame(x1 = runif(n, -3, 3),
  id = as.factor(rep(names(MunichBnd), length.out = n)))
dat$sp <- with(dat, sort(runif(N, -2, 2), decreasing = TRUE)[id])
     
## response
dat$y <- with(dat, 1.5 + sin(x1) + sp + rnorm(n, sd = 0.6))

## sort data according id
dat <- dat[order(x2int(dat$id)),]

## estimate models with
## BayesX MCMC and mgcv REML
xt <- list(polys = MunichBnd)
b1 <- bayesx(y ~ s(x1, bs = "ps", k = 10) + 
  s(id, bs = "mrf", xt = xt), 
  method = "MCMC", data = dat, dir.rm = FALSE)
b2 <- gam(y ~ s(x1, bs = "ps", k = 10) + 
  s(id, bs = "mrf", xt = xt), 
  method = "REML", data = dat)

plot(b1, term = "s(x1)", 
  c.select = c("x1", "Mean", "2.5%", "97.5%"),
  fill.select = c(0, 0, 1, 1))

p2 <- predict(b2, terms = 2)
par(mfrow = c(1, 3))
plot(b1, term = "s(id)", map = MunichBnd, symmetric = FALSE, range= c(-3, 3))
plotmap(MunichBnd, x = p2, id = names(p2), symmetric = FALSE, range= c(-3, 3))
plotmap(MunichBnd, x = unique(cbind(x2int(dat$id), dat$sp)), 
  main = "Truth",  range= c(-3, 3))




## zambia data
# dat <- read.table("http://www.stat.uni-muenchen.de/~bayesx/tutorials/zambia.raw", header = TRUE)
ZambiaNutrition <- list()
ZambiaNutrition$stunting <- dat$hazstd
ZambiaNutrition$bmi <- dat$bmi
ZambiaNutrition$agechild <- dat$agc
ZambiaNutrition$district <- dat$district
memployment <- rep("working", nrow(dat))
memployment[dat$rcw == -1] <- "not working"
ZambiaNutrition$memployment <- as.factor(memployment)
education <- rep("no/incomplete", nrow(dat))
education[dat$edu1 == 1] <- "minimum primary"
education[dat$edu2 == 1] <- "minimum secondary"
ZambiaNutrition$education <- as.factor(education)
urban <- rep("no", nrow(dat))
urban[dat$tpr == 1] <- "yes"
ZambiaNutrition$urban <- as.factor(urban)
gender <- rep("male", nrow(dat))
gender[dat$sex == -1] <- "female"
ZambiaNutrition$gender <- gender
ZambiaNutrition <- as.data.frame(ZambiaNutrition)
ZambiaNutrition <- ZambiaNutrition[order(ZambiaNutrition$district),]
summary(ZambiaNutrition)
## save(ZambiaNutrition, file = "/home/nikolaus/svn/bayesr/pkg/R2BayesX/data/ZambiaNutrition.rda")


## forest health data
# dat <- as.matrix(read.table("http://www.stat.uni-muenchen.de/~kneib/regressionsbuch/download/buche.raw", header = TRUE))
dat[dat == "."] <- NA
storage.mode(dat) <- "numeric"
dat <- as.data.frame(dat)
ForestHealth <- list()
ForestHealth$id <- dat$id 
ForestHealth$year <- dat$jahr 
buche3 <- rep(0, nrow(dat))
buche3[dat$buche >= 12.5] <- 1 
buche3[dat$buche >= 50] <- 2
ForestHealth$defoliation <- buche3 # dat$buche
ForestHealth$x <- dat$x
ForestHealth$y <- dat$y
ForestHealth$age <- dat$alter 
ForestHealth$canopy <- dat$schirm 
ForestHealth$inclination <- dat$hang 
ForestHealth$elevation <- dat$hoehe
ForestHealth$soil <- dat$grund 
ForestHealth$ph <- dat$ph 
ForestHealth$moisture <- as.factor(dat$frische) 
ForestHealth$alkali <- as.factor(dat$alkali) 
H <- rep(0, nrow(dat))
H[dat$humus == 1] <- 1
H[dat$humus == 2] <- 2
H[dat$humus == 3] <- 3
H[dat$humus > 3] <- 4
ForestHealth$humus <- as.factor(H)
ForestHealth$stand <- as.factor(dat$artkat)
levels(ForestHealth$stand) <- c("mixed", "deciduous")
ForestHealth$fertilized <- as.factor(dat$dueng)
levels(ForestHealth$fertilized) <- c("no", "yes")
ForestHealth <- na.omit(as.data.frame(ForestHealth))
## save(ForestHealth, file = "/home/nikolaus/svn/bayesr/pkg/R2BayesX/data/ForestHealth.rda")


## large data example
library("R2BayesX")
n <- 100000
file <- paste(tempdir(), "/ldata.raw", sep = "")
write.table(matrix(c("x", "y"), nrow = 1), file = file, 
  quote = FALSE, row.names = FALSE, col.names = FALSE)
for(i in 1:10) {
  dat <- data.frame(x = round(runif(n, -3, 3), 2L))
  dat$y <- with(dat, sin(x) + rnorm(n, sd = 2))
  write.table(dat, file = file, append = TRUE, 
    quote = FALSE, row.names = FALSE, col.names = FALSE)
}
b <- bayesx(y ~ s(x, bs = "ps"), data = file, 
  iter = 3000, burnin = 500, step = 2, predict = FALSE)


## create a package skeleton
dir <- "/home/c403129/svn/bayesr/pkg/R2BayesX/R"
path <- tempdir()
invisible(sapply(paste(dir, "/", list.files(dir), sep = ""), source))
Rfiles <- ls(all.names = TRUE)
Rfiles <- Rfiles[!Rfiles %in% c("dir", "path", "Rfiles")]
package.skeleton(name = "r2bayesx", Rfiles, path = path)







