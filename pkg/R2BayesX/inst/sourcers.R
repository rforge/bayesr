dir <- "/home/c403129/svn/bayesr/pkg/R2BayesX/R"
## dir <- "J:/c403/stat/R2BayesX/R"
invisible(sapply(paste(dir, "/", list.files(dir), sep = ""), source))
plot(b, which = "coef-samples", all.acf = TRUE)


script <- getscript(b2)




sc <- getscript(b2, file = "~/tmp/getscript.R")

b <- read.bayesx.output("~/tmp/fit4")





fit_z_new1<-bayesx(stunt ~ -1+c_sex + residence0 + residence1 + residence2 + precare +
  bornhome + fhh +
  s(age_c,k=22,bs="ps",xt=list(lambda=1000,binning=50,centermethod="meansimple")) +
  s(age_c,k=22,bs="ps",by=c_sex,xt=list(lambda=1000,binning=50,centermethod="meansimple")) +
  s(bmi_mo_c,k=22,bs="ps",xt=list(lambda=1000,binning=50,centermethod="meansimple")) +
  s(biage_c,k=22,bs="ps",xt=list(lambda=1000,binning=50,centermethod="meansimple")) +
  s(vac_numb_c,k=8,bs="ps",xt=list(lambda=1000,binning=50,centermethod="meansimple")) +
  s(bi_pre_c,k=22,bs="ps",xt=list(lambda=1000,binning=50,centermethod="meansimple")) +
  s(ai_distdiff,k=22,bs="ps",xt=list(lambda=1000,binning=50,centermethod="meansimple")) +
  s(educm_y_distdiff,k=12,bs="ps",xt=list(lambda=1000,binning=50,centermethod="meansimple")) +
  s(hhs_c,k=12,bs="ps",xt=list(lambda=1000,binning=50,centermethod="meansimple")) +
  s(b_order_c,k=8,m=c(3,2),bs="ps",xt=list(lambda=1000,binning=50,centermethod="meansimple")) +
  r(distH, ~ -1 +
  s(dist_ai_c,k=22,bs="ps",xt=list(lambda=1000,binning=50,centermethod="meansimple"))+
  s(dist_eduyears_c,k=12,bs="ps",xt=list(lambda=1000,binning=50,centermethod="meansimple")) +
  r(region,~s(gdp_c,k=12,bs="ps",xt=list(lambda=1000,binning=50,centermethod="meansimple")), 
  data=dregion), data=ddist), data=d , verbose=FALSE, method="HMCMC", iter=3000, step=1, burnin=2000, 
  family="binomial_logit", dir.rm=FALSE)





## now estimate a random effects
## model with a stage 2 covariate
## effect
set.seed(333)
n <- 1000

## 1st stage
N <- 500
dat1 <- data.frame(id = sort(rep(1:N, n/N)), 
  x1 = runif(n, -3, 3))
dat1$re <- with(dat1, rnorm(N, sd = 0.6)[id])

## 2nd stage
N2 <- 100
dat2 <- data.frame(id = unique(dat1$id), 
  x2 = runif(N, -3, 3), id2 = sort(rep(1:N2, N/N2)))
dat2$re2 <- with(dat2, rnorm(N2, sd = 0.6)[id2])

## 3rd stage
dat3 <- data.frame(id2 = unique(dat2$id2), 
  x3 = runif(N2, 0, 5))

## response
dat1$y <- with(dat1, 1.5 + sin(x1) + re + cos(dat2$x2)[id] + 
  dat2$re2[id] + dat3$x3[dat2$id2][id] + rnorm(n, sd = 0.6))

## estimate hierarchical model
## with the intercept in the 
## 2nd stage
system.time(b1 <- bayesx(y ~ -1 + s(x1, bs = "ps") + 
  r(id, ~ -1 + s(x2 , bs = "ps") + r(id2, ~ 1 + s(x3, bs = "ps"), data = dat3), data = dat2), 
  method = "HMCMC", data = dat1, iter = 3000, burnin = 1000, outfile = "~/tmp/bayesxH"))



B <- read.bayesx.output("/tmp/RtmpePekC5/bayesx")






m1 <- bayesx(y ~ -1 + s(x1) + r(id1 ~ -1 + s(r1) + r(id2 ~ -1 + s(r2) + r(id3 ~1 + s(r3)))), dir.rm = FALSE)





     system.time(b1 <- bayesx(y ~ -1 + s(x1, bs = "ps") + 
       r(id ~ 1 + s(x2 , bs = "ps"), data = dat2), 
       method = "HMCMC", data = dat1, iter = 3000, burnin = 1000, step = 2, dir.rm = FALSE))



B <- read.bayesx.output("/tmp/RtmpQ9qdbE/bayesx")






p <- parse.bayesx.input(y ~ -1 + s(x1, bs = "ps") + 
  r(id, ~ 1 + s(x2 , bs = "ps") + r(id3, ~ s(x2)), data = dat2, seed = 333), 
  method = "MCMC", data = dat1, seed = 111)
write.bayesx.input(p)




BIC(c(fm1, fm2))

b <-c(fm1, fm2)



ff <- stunting ~ memployment + education + urban + gender
m <- as.data.frame(model.matrix(ff, data = ZambiaNutrition, contrasts.arg = list(memployment = "contr.sum", 
  education = "contr.sum", urban = "contr.sum", gender = "contr.sum")))
cbind(m["memployment1"], ZambiaNutrition["memployment"])

b <- read.bayesx.output("/tmp/RtmpvqFbDK/bayesx")






par(xaxs = "i", yaxs = "i")
plot(0, 1, type = "n", xlim = c(0, 1), ylim = c(0, 1))
rect(0:98/99, 0, 1:99/99, 1, border = NA, 
  col = make_pal(col = diverge_hcl, ncol = 99, range = c(0, 1))$colors)



     
## generate some data
set.seed(111)
n <- 200
     
## regressor
dat <- data.frame(x = runif(n, -3, 3))
     
## response
dat$y <- with(dat, 1.5 + sin(x) + rnorm(n, sd = 0.6))
     
## estimate models with
## bayesx REML and MCMC
b <- bayesx(y ~ s(x, bs = "ps"), 
  method = "HMCMC", data = dat)

bayesx_logfile(b)








summary(b)


b <- bayesx(mstatus ~ s(age, by = myvar, bs = "ps", k = 20), method = "REML",
  family = "multinomial", data = nzmarital, reference = 2)





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
set.seed(123)
n <- 100000
file <- paste(tempdir(), "/data.raw", sep = "")
write.table(matrix(c("x", "y"), nrow = 1), file = file, 
  quote = FALSE, row.names = FALSE, col.names = FALSE)
for(i in 1:50) {
  dat <- data.frame(x = round(runif(n, -3, 3), 2L))
  dat$y <- with(dat, sin(x) + rnorm(n, sd = 2))
  write.table(dat, file = file, append = TRUE, 
    quote = FALSE, row.names = FALSE, col.names = FALSE)
}
b <- bayesx(y ~ s(x, bs = "ps"), data = file, 
  iter = 3000, burnin = 500, step = 2, predict = FALSE)
save(b, file = "/home/c403129/svn/bayesr/pkg/R2BayesX/inst/JSS/large-data.rda")


dat <- read.table(file, header = TRUE)
b <- gam(y ~ s(x, bs = "ps"), data = dat)

## create a package skeleton
dir <- "/home/c403129/svn/bayesr/pkg/R2BayesX/R"
path <- tempdir()
invisible(sapply(paste(dir, "/", list.files(dir), sep = ""), source))
Rfiles <- ls(all.names = TRUE)
Rfiles <- Rfiles[!Rfiles %in% c("dir", "path", "Rfiles")]
package.skeleton(name = "r2bayesx", Rfiles, path = path)


## multinomial example
library("R2BayesX")
data("nzmarital", package = "VGAM")
b <- bayesx(mstatus ~ s(age, bs = "ps"), method = "REML",
  family = "multinomial", data = nzmarital, reference = "Married/Partnered")
fb <- fitted(b)[,4:6]
cnfb <- c(gsub("mu:", "", colnames(fb)), "Married/Partnered")
fb <- cbind(fb, 1 - rowSums(fb))
mycol <- c("red", "darkgreen", "blue")
ooo <- with(nzmarital, order(age))
with(nzmarital, matplot(age[ooo], fb[ooo,],
  type = "l", las = 1, lwd = 2, ylim = 0:1, ylab = "Fitted probabilities",
  xlab = "Age", col = c(mycol[1], "black", mycol[-1])))
legend(x = 52.5, y = 0.62, col = c(mycol[1], "black", mycol[-1]),
  lty = 1:4, legend = cnfb, lwd = 2)


fit.ms <- vgam(mstatus ~ s(age, df = 10), multinomial(refLevel = 2), data = nzmarital)
par(mfrow = c(3, 1))
plotvgam(fit.ms, se = TRUE, which.term = 1)
plotvgam(fit.ms, se = TRUE, which.term = 2)
plotvgam(fit.ms, se = TRUE, which.term = 3)
mu <- fitted(fit.ms)


## stepwise example
## an example of automatic model selection via null space penalization
set.seed(3); n <- 500
dat <- gamSim(1, n = n, scale = 1.2) ## simulate data
dat$x4 <- runif(n, 0, 1); dat$x5 <- runif(n, 0, 1) ## spurious

b <- gam(y ~ s(x0, bs = "ps") + s(x1, bs = "ps") + s(x2, bs = "ps") + 
  s(x3, bs = "ps") + s(x4, bs = "ps") + s(x5, bs = "ps"),
  data = dat, select = TRUE, method = "REML")

summary(b)
plot(b, pages = 1)

b2 <- bayesx(y ~ s(x0, bs = "ps") + s(x1, bs = "ps") + s(x2, bs = "ps") + 

  s(x3, bs = "ps") + s(x4, bs = "ps") + s(x5, bs = "ps"),
  data = dat, family = "gaussian", method = "STEP", dir.rm = FALSE)

b3 <- bayesx(y ~ s(x0, bs = "ps", k = 20) + s(x1, bs = "ps", k = 20) + s(x2, bs = "ps", k = 20),
  data = dat, family = "gaussian", method = "REML")


## a dummy example
## more examples
set.seed(111)
n <- 500
     
## regressors
dat <- data.frame(x = runif(n, -3, 3), z = runif(n, -3, 3),
  w = runif(n, 0, 6), fac = factor(rep(1:2, n/2)))
     
## response
dat$y <- with(dat, 1.5 + sin(x) + cos(z) * sin(w) +
  c(2.67, 3.33)[fac] + rnorm(n, sd = 0.6))

b1 <- bayesx(y ~ s(x, bs = "ps") + s(z, w, bs = "te") + fac,
  data = dat, method = "MCMC")

b2 <- bayesx(y ~ s(x, bs = "ps") + s(z, w, bs = "te") + fac,
  data = dat, method = "MCMC", contrasts = list(fac = contr.sum))







