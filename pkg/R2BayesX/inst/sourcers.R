dir <- "/home/nikolaus/svn/bayesr/pkg/R2BayesX/R"
## dir <- "J:/c403/stat/R2BayesX/R"
invisible(sapply(paste(dir, "/", list.files(dir), sep = ""), source))
plot(b1, which = "scale-samples")




    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")





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
ForestHealth$id <- dat$id 
ForestHealth$year <- dat$jahr 
ForestHealth$beech <- dat$buche
ForestHealth$x <- dat$x
ForestHealth$y <- dat$y
ForestHealth$age <- dat$alter 
ForestHealth$canopy <- dat$schirm 
ForestHealth$inclination <- dat$hang 
ForestHealth$elevation <- dat$hoehe
ForestHealth$soil <- dat$grund 
ForestHealth$ph <- dat$ph 
ForestHealth$moisture <- as.factor(dat$frische) 
ForestHealth$saturation <- as.factor(dat$alkali) 
ForestHealth$humus <- dat$humus
ForestHealth$stand <- as.factor(dat$artkat)
levels(ForestHealth$stand) <- c("deciduous", "mixed")
ForestHealth$fertilized <- as.factor(dat$dueng)
levels(ForestHealth$fertilized) <- c("yes", "no")
ForestHealth <- as.data.frame(ForestHealth)
## save(ForestHealth, file = "/home/c403129/svn/bayesr/pkg/R2BayesX/data/ForestHealth.rda")













