dir1 <- "/home/c403129/svn/bayesr/pkg/R2BayesX/R"
## dir1 <- "J:/c403/stat/R2BayesX/R"
files <- list.files(dir1)
for(i in 1:length(files))
  source(paste(dir1, "/", files[i], sep = "")) 
b <- read.bayesx.output("/tmp/RtmpjoTH9E/bayesx")






colorlegend(side.legend = 2, pos = "topleft", xpd = TRUE)






b1 <- bayesx(y ~ s(x, bs = "ps"), method = "REML", data = dat, dir.rm = FALSE)

 

data("ZambiaNutrition")
data("ZambiaBnd")
b <- bayesx(stunting ~ memployment + education + urban + gender + 
  s(bmi, bs = "ps") + s(agechild, bs = "ps") +
  s(district, bs = "mrf", xt = list(map = ZambiaBnd)) + r(district),
  iter = 1200, burnin = 200, step = 1, data = ZambiaNutrition)









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
## save(ZambiaNutrition, file = "/home/c403129/svn/R2BayesX/data/ZambiaNutrition.rda")


b <- bayesx(y ~ s(x, bs = "ps", xt = list(contourprob = 4)))


plot(b, which = 1)

   




txt <- "
olkdlklklkdlk: is here
      Default: not here!
"



k <- 100
zk <- seq(min(dat$z), max(dat$z), length = k)
wk <- seq(min(dat$w), max(dat$w), length = k)
knots <- cbind(zk, wk)

b <- bayesx(y ~ s(x, bs = "ps") + s(z, w, bs = "kr", xt = list(knotdata = knots)) + fac, 
  method = "REML", data = dat)


plot(b, resid = F, image = T)
















barplot(rep(1, 9), col = diverge_hcl(9), axes = FALSE, border = NA)



plot(b, term = 2, map = munich)







plot(b, term = "s(id)", map = munich)






b1 <- bayesx(y ~ s(x, bs = "ps") + te(z, w) + fac,
        data = dat, method = "MCMC")

plot(b, term = "s(id)")

par(mfrow = c(1, 2))
plot(b1, term = 2, map = munich, main = "Estimate", values = TRUE, range = c(-2.5, 2.5))
plotmap(munich, x = unique(cbind(dat$id, dat$sp)), main = "Truth", values = TRUE, range = c(-2.5, 2.5))





plot(b, term = "s(id)", map = munich, range = c(-1, 1))



     id       x       y   Mean
1     1 4468800 5333390  2.554
2     2 4468530 5333000  2.554
3     3 4468030 5333150  2.479
4     4 4468180 5333580  2.560
5     5 4469300 5333370  2.356
6     6 4469590 5334820  1.554
7     7 4468660 5332860  2.517
8     8 4469030 5332610  2.404
9     9 4468060 5332160  1.952
10   10 4467560 5331480  1.516
11   11 4467600 5332410  2.063
12   12 4466900 5331850  1.929
13   13 4467420 5333270  2.240
14   14 4466420 5332860  1.730
15   15 4467950 5334250  2.383
16   16 4467390 5334360  2.229
17   17 4466850 5334820  1.916
18   18 4466390 5334220  1.868
19   19 4467710 5335280  2.036
20   20 4468190 5335300  2.055
21   21 4468690 5334680  2.133
22   22 4469170 5334690  1.889
23   23 4467230 5335200  1.917
24   24 4468180 5335670  1.969
25   25 4468460 5337190  1.649
26   26 4467260 5336020  1.705
27   27 4470000 5333560  1.821
28   28 4470760 5333440  1.439
29   29 4470920 5332860  1.571
30   30 4470370 5332070  1.610
31   31 4469520 5331840  1.635
32   32 4468850 5331780  1.466
33   33 4465660 5331040  1.608
34   34 4466810 5330790  1.257
35   35 4465270 5330320  1.218
36   36 4464890 5332710  0.397
37   37 4463690 5330470  1.060
38   38 4465590 5333350  1.185
39   39 4466310 5332850  1.680
40   40 4464990 5335440  1.098
41   41 4462590 5335920  1.229
42   42 4466460 5336590  1.310
43   43 4465940 5334570  1.551
44   44 4466550 5335260  1.581
45   45 4465380 5336190  1.074
46   46 4464940 5338540  0.954
47   47 4462600 5338050  0.729
48   48 4468840 5341870  0.270
49   49 4466880 5338370  1.342
50   50 4468670 5338430  1.018
51   51 4470740 5341370  0.927
52   52 4473260 5341390  1.284
53   53 4471020 5337080  0.796
54   54 4469290 5336760  1.036
55   55 4470070 5336210  0.599
56   56 4469320 5335260  1.497
57   57 4470330 5335570  0.712
58   58 4469390 5338240  0.596
59   59 4472500 5337490  1.079
60   60 4475730 5337170  0.665
61   61 4471290 5335620  0.687
62   62 4472690 5334010  0.999
63   63 4475050 5334570  0.921
64   64 4471910 5334220  0.923
65   65 4470430 5333780  1.455
66   66 4472420 5332010  0.400
67   67 4475620 5334050  0.455
68   68 4478360 5332560  0.172
69   69 4474460 5330650  0.239
70   70 4477330 5329530 -0.257
71   71 4471550 5331070  0.070
72   72 4470430 5330010  0.586
73   73 4472650 5327770 -0.400
74   74 4474060 5329670  0.033
75   75 4475770 5327320 -0.437
76   76 4469360 5330490  0.138
77   77 4470230 5328070 -0.118
78   78 4467980 5331070  0.650
79   79 4466900 5329340  0.027
80   80 4468500 5329910 -0.274
81   81 4469060 5328760  0.207
82   82 4466560 5327000 -0.289
83   83 4465350 5326770 -0.303
84   84 4466250 5328670 -0.504
85   85 4464760 5328860 -0.527
86   86 4462340 5327070 -0.547
87   87 4461150 5327830 -0.606
88   88 4463550 5326360 -0.687
89   89 4461060 5332230 -1.125
90   90 4460770 5331400 -1.006
91   91 4460700 5330010 -0.877
92   92 4460890 5334970 -0.630
93   93 4461590 5333720 -1.351
94   94 4459340 5333000 -0.673
95   95 4460910 5336680 -0.866
96   96 4454420 5337380 -1.284
97   97 4454770 5334310 -0.486
98   98 4455800 5338800 -1.573
99   99 4461300 5340190 -0.987
100 100 4458900 5339450 -0.825
101 101 4463750 5344090 -1.018
102 102 4466830 5340790 -1.137
103 103 4462360 5341020 -1.164
104 104 4464880 5340100 -1.284
105 105 4464160 5333060 -0.693
106 106 4462470 5333350 -1.566




make.legend(x = x, at = c(0, 0.5, 1))


plotmap(germany)

plotmap(munich, x = x1[[1]], which = "90%")









b <- read.bayesx.output("/tmp/RtmptYKZo9/bayesx", "bayesx.estim")



m <- read.bayesx.output(file)


dir <- "/home/c403129/Desktop/bayesx.test.files"
# dir <- "/home/nikolaus/svn/bayesx.test.files"

file <- paste(dir,"/bayesx_by2",sep="")dir1 <- "/home/c403129/svn/R2BayesX/R"
files <- list.files(dir1)
for(i in 1:length(files))
  source(paste(dir1,"/",files[i],sep=""))
b <- bayesx(y ~ s(x1, bs = "ps") + 
  r(id, ~ 1 + s(id, bs = "mrf", xt = list(map = munich)), data = dat2), 
  method = "MCMC", data = dat1)





m <- read.bayesx.output(file, "bayesx.estim")

file <- paste(dir,"/bayesx_hier",sep="")
m <- read.bayesx.output(file, "bayesx.estim_hlevel1_MAIN_REGRESSION")

file <- paste(dir,"/bayesx_mcmc",sep="")
m <- read.bayesx.output(file, "bayesx.estim")

file <- paste(dir,"/bayesx_reml",sep="")
m <- read.bayesx.output(file, "bayesx.estim")

file <- paste(dir,"/bayesx_step",sep="")
m <- read.bayesx.output(file, "bayesx.estim")

file <- paste(dir,"/bayesx_geo",sep="")
m <- read.bayesx.output(file, "bayesx.estim")

file <- paste(dir,"/bayesx_by",sep="")
m <- read.bayesx.output(file, "bayesx.estim")

file <- paste(dir,"/bayesx_hier_by",sep="")
m <- read.bayesx.output(file, "bayesx.estim_hlevel3_RANDOM_EFFECTS")






## Not run:
     
## generate some data
set.seed(111)
n <- 500
## regressors
ex1 <- data.frame(x_one = runif(n, -3, 3), z_three = runif(n, -3, 3),
  w_two = runif(n, 0, 6), fac = factor(rep(1:10, n/10)))
## response
ex1$y <- with(ex1, 10 + sin(x_one) + cos(z_three) * sin(w_two) +
  c(2.67, 5, 6, 3, 4, 2, 6, 7, 9, 7.5)[fac] + rnorm(n, sd = 0.6))
     
## estimate models with
## bayesx MCMC and REML
## and compare with
## gam() GCV
b1 <- bayesx(y ~ s(x_one, bs = "ps") + te(z_three, w_two, bs = "ps") + fac,
  data = ex1, method = "MCMC", dir.rm = FALSE, outfile = "/home/c403129/Desktop/bayesx.test.files/bayesx_mcmc", 
  iter = 12000, burnin = 2000, step = 10)

b2 <- bayesx(y ~ s(x_one, bs = "ps") + te(z_three, w_two, bs = "ps") + fac,
  data = ex1, method = "REML", dir.rm = FALSE, outfile = "/home/c403129/Desktop/bayesx.test.files/bayesx_reml")

b3 <- bayesx(y ~ s(x_one, bs = "ps") + te(z_three, w_two, bs = "ps") + fac,
  data = ex1, method = "STEP", dir.rm = FALSE, outfile = "/home/c403129/Desktop/bayesx.test.files/bayesx_step")

b4 <- bayesx(y ~ s(x_one, bs = "ps") + s(z_three, by = w_two) + fac,
  data = ex1, method = "MCMC", dir.rm = FALSE, outfile = "/home/c403129/Desktop/bayesx.test.files/bayesx_by")

b5 <- bayesx(y ~ s(x_one, bs = "ps") + s(z_three, by = fac),
  data = ex1, method = "MCMC", dir.rm = FALSE, outfile = "/home/c403129/Desktop/bayesx.test.files/bayesx_by2")







s <- search.bayesx.models("/tmp/RtmpLGdcTk/bayesx_hier")

b1 <- bayesx(y ~ s(x_one, bs = "ps") + s(z_three,by=w_two) + fac, method = "REML", data = ex1, dir.rm = FALSE)

p <- parse.bayesx.input(y~-1+s(x1)+r(id1,~-1+s(x2)+r(id2,~1+s(x3))))
b <- write.bayesx.input(p)

s3 <- bayesx(y~-1+s(x1)+r(id1,~-1+s(x2)+r(id2,~1+s(x3))), dir.rm=FALSE)









b1 <- bayesx(y ~ s(x, bs = "ps") + te(z,w) + fac, method = "MCMC", data = ex1)



m <- model.frame(formula = y ~ fac + z + w + x + z + w, data = ex1, weights = w, offset = z, drop.unused.levels = TRUE)




b1 <- bayesx(y ~ s(x, bs = "ps") + te(z,w) + fac, method = "REML", data = ex1, dir.rm = FALSE)
b1 <- bayesx(y ~ s(x, bs = "ps") + te(z,w) + fac + z + w, method = "MCMC", data = ex1)
b1 <- bayesx(y ~ s(x, bs = "ps") + fac + s(z) + s(w), method = "MCMC", data = ex1)
b1 <- bayesx(y ~ s(x, bs = "ps") + te(z,w) + fac, method = "REML", data = ex1, dir.rm=F)
b1 <- bayesx(y~s(x1,bs="ps")+s(id,bs="mrf",xt=list(map=munich)))
b1 <- read.bayesx.output("/tmp/Rtmp8lqSkV/bayesx", model.name = "bayesx.estim")
b1 <- bayesx(y ~ s(x, bs = "ps") + te(z,w) + fac, method = "REML", data = ex1)
b1 <- read.bayesx.output("/tmp/RtmpQHsy8j/bayesx", model.name = "bayesx.estim")







model_f_x_pspline.res
model_MAIN_REGRESSION_nonlinear_pspline_effect_of_x1.res


b1 <- bayesx(y~-1+s(x1,bs="ps")+r(id,~1+s(x2,bs="ps")))




library("R2BayesX")


### BAYESX INSTALLING
check.install.bayesx()
install.bayesx()


### BAYESX EXAMPLES
# set BayesX bin path
options(bayesx.bin="/home/c403129/bin/BayesX")

# generate some data
set.seed(111)
n<-500
x<-runif(n,-3,3); z<-runif(n,-3,3); w<-runif(n,0,6)
fac<-rep(1:10,n/10)

# generate some nonlinear function 1d and 2d
f1<-function(x){sin(x)}
f2<-function(x,z){cos(x)*sin(z)}

# some values for the factor covariate
betas<-c(2.67,5,6,3,4,2,6,7,9,7.5)[fac]
fac<-as.factor(fac)

# generate response
y<-10+f1(x)+f2(z,w)+betas+rnorm(n,sd=0.6)

# estimate models with
# bayesx MCMC and REML
# and compare with 
# gam() GCV
b1<-bayesx(y~s(x,bs="ps")+te(z,w,bs="ps")+fac)
b2<-bayesx(y~s(x,bs="ps")+te(z,w,bs="ps")+fac,method="REML")
b3<-gam(y~s(x,bs="nu")+te(z,w,bs="nu")+fac)

# summary statistics
summary(b1)
summary(b2)
summary(b3)

# plot the effects
par(mfrow=c(3,2),mar=c(2,2,2,0))
plot(b1,1); plot(b1,2,ticktype="simple")
plot(b2,1); plot(b2,2,ticktype="simple")
plot(b3,select=1); vis.gam(b3,c("z","w"),theta=40,phi=40)

# objects in b
names(b1)


### BAYESX OUTPUT HANDLING FUNCTIONS
# read in all model output from folder
models<-read.bayesx.output("/home/nikolaus/svn/r-test/pspline/output")
print(models)
summary(models)
plot(models)

# read in only specified model from folder
b<-read.bayesx.output("/home/nikolaus/svn/r-test/pspline/output","step_gauss_gam1")
summary(b)
plot(b,all=TRUE)

# 2nd example
# first read in map of tanzania
tanzania <- read.bnd("/home/nikolaus/svn/R2BayesX/inst/examples/ex02/tanzania.bnd")

b<-read.bayesx.output("/home/nikolaus/svn/R2BayesX/inst/examples/ex02","res")
summary(b)
plot(b,all=TRUE,map=tanzania)

# 3rd example
germany <- read.bnd("/home/nikolaus/svn/R2BayesX/inst/examples/ex03/germany.bnd")
b<-read.bayesx.output("/home/nikolaus/svn/R2BayesX/inst/examples/ex03","spatial")
summary(b)
plot(b,map=germany,lpos=c(0.7,0.5),pal="heat_hcl")

# last one
b<-read.bayesx.output("/home/nikolaus/svn/R2BayesX/inst/examples/ex03","surface")
summary(b)
plot(b)
# need to check swaping???
plot(b,swap=T)






## hierarchical
# setup a model with 2 random effect stages
set.seed(333)
n<-2000
N<-1000
x1<-runif(n,-3,3)
id1<-as.factor(round(runif(n,1,N)))
nr<-nlevels(id1)
x2<-runif(nr,-3,3)

NN<-500
id2<-as.factor(round(runif(nr,1,NN)))
nr2<-nlevels(id2)
x3<-runif(nr2,-2,4)
x4<-runif(nr2,-3,3)

# some functions for the terms
fr2<-cos(x3) + sin(x4)
fr2<-fr2-mean(fr2)
re2<-rnorm(nr2,sd=0.6)
re2<-re2-mean(re2)

fr1<-sin(x2)
fr1<-fr1-mean(fr2)
re1<-rnorm(nr,sd=0.6)
re1<-re1-mean(re1)
random<-c(fr1+re1+c(fr2+re2)[id2])[id1]
y<-cos(x1)+random+rnorm(n,sd=0.6)

d1 <- data.frame(y=y,x1=x1,id1=id1)
d2 <- data.frame(id1=unique(id1),x2=x2,id2=id2)
d3 <- data.frame(id2=unique(id2),x3=x3,x4=x4)


m <- bayesx(y~-1+s(x1)+r(id1,~-1+s(x2)+r(id2,~1+s(x3) + s(x4), data = d3), data = d2), data = d1, 
  outfile = "/home/c403129/Desktop/bayesx.test.files/bayesx_hier")

m <- bayesx(y~s(x1,bs="ps")+r(id1), data = d1, dir.rm = FALSE)


p <- parse.bayesx.input(y~-1+s(x1)+r(id1,~-1+s(x2)+r(id2,~1+s(x3) + s(x4), data = d3), data = d2), data = d1, 
  outfile = "/home/c403129/Desktop/bayesx.test.files/bayesx_hier")

write.bayesx.input(p)

