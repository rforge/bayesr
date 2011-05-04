dir <- "/home/nikolaus/svn/bayesr/pkg/R2BayesX/R"
## dir <- "J:/c403/stat/R2BayesX/R"
invisible(sapply(paste(dir, "/", list.files(dir), sep = ""), source))
data("ForestHealth", package = "R2BayesX")
data("BeechBnd")
fm2 <- bayesx(defoliation ~  stand + fertilized + 
  humus + moisture + alkali + ph + soil + 
  s(age, bs = "ps") + s(inclination, bs = "ps") + 
  s(canopy, bs = "ps") + s(year, bs = "ps") + 
  s(elevation, bs = "ps") + 
  s(x, y, bs = "te"),
  family = "cumlogit", method = "REML", data = ForestHealth)





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













