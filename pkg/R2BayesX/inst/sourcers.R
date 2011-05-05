dir <- "/home/c403129/svn/bayesr/pkg/R2BayesX/R"
## dir <- "J:/c403/stat/R2BayesX/R"
invisible(sapply(paste(dir, "/", list.files(dir), sep = ""), source))

set.seed(333)
     
## simulate some geographical data
data("MunichBnd")
N <- length(MunichBnd); names(MunichBnd) <- 1:N
n <- N*5
     
## regressors
dat <- data.frame(id = rep(1:N, n/N), x1 = runif(n, -3, 3))
dat$sp <- with(dat, sort(runif(N, -2, 2), decreasing = TRUE)[id])
     
## response
dat$y <- with(dat, 1.5 + sin(x1) + sp + rnorm(n, sd = 0.6))
     
## estimate models with
## BayesX MCMC and REML
xt <- list(polys = MunichBnd)
b1 <- bayesx(y ~ s(x1, bs = "ps", k = 10) + 
  s(id, bs = "mrf", xt = xt), 
  method = "MCMC", data = dat)
b2 <- gam(y ~ s(x1, bs = "ps", k = 10) + 
  s(id, bs = "mrf", xt = xt), 
  method = "REML", data = dat)

plot(b1, c.select = c("x1", "Mean", "2.5%", "97.5%"),
  fill.select = c(0, 0, 1, 1))

par(mfrow = c(2, 1))
plot(b1, term = "s(id)", map = MunichBnd, col = heat.colors)
plot(b2, scheme = "heat", select = 2)
plotmap(MunichBnd, x = fitted(b2, type = "terms"), dat$id)




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













