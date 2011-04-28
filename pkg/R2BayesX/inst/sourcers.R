dir <- "/home/c403129/svn/bayesr/pkg/R2BayesX/R"
## dir <- "J:/c403/stat/R2BayesX/R"
invisible(sapply(paste(dir, "/", list.files(dir), sep = ""), source))
b1 <- bayesx(y ~ s(x, bs = "ps"), data = dat)









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
# dat <- read.table("http://www.stat.uni-muenchen.de/~kneib/regressionsbuch/download/buche.raw", header = TRUE)
ForestHealth <- list()
ForestHealth$year <- dat$jahr 
ForestHealth$x <- dat$x
ForestHealth$y <- dat$y
ForestHealth$slope <- dat$hang 
ForestHealth$altitude <- dat$hoehe
ForestHealth$base <- dat$grund 
ForestHealth$fertilization <- dat$dueng
ForestHealth$age <- dat$alter 
ForestHealth$shield <- dat$schirm 
ForestHealth$artkat <- dat$artkat
ForestHealth$humus <- dat$humus 
ForestHealth$beech <- dat$buche
ForestHealth$id <- dat$id 
ForestHealth$alkali <- dat$alkali 
ForestHealth$ph <- dat$ph 
ForestHealth$viridity <- dat$frische 
ForestHealth <- as.data.frame(ForestHealth)
## save(ForestHealth, file = "/home/nikolaus/svn/bayesr/pkg/R2BayesX/data/ForestHealth.rda")













