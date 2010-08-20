\name{gibbs}
\alias{gibbs}

\title{Bayesian Hierarchical Generalized Structured Additive Regression Models}

\description{This routine fits a bayesian hierarchical generalized structured additive regression model to the data. Nonlinear smooth functions of covariates 
can be modeled with b/p-splines, tensor splines, univariate and bivariate kriging methods. Spatial heterogeneity can be modeled with markov random fields and/or random effects constructors

Within the model framework it is possible to model random effects in terms of any covariate available and function as described above. Random effects may 
be modeled hierarchical or non hierarchical. It is also possible to model various varying effects on different levels through any of the nonlinear smooth functions

Hyperparameters for the smooth functions are estimated automatically within the bayesian framework
}
 
\usage{
gibbs(formula, data, weights, family = "gaussian", iter = 1200, 
      burnin = 200, thinning = 1, hyperasigma = 1e-04, hyperbsigma = 1e-04, 
      hyperatau = 1e-04, hyperbtau = 1e-04, trace = TRUE, dots = FALSE, 
      speed = TRUE, ...) 
}

\arguments{
  \item{formula}{a \code{\link{formula}} type object, where the response variable is on the left side of an ~ operator, 
                 every term is seperated by the + operator on the right hand side. Nonlinear smooth functions, random and multilevel effects,
		             which should enter the model are simply integrated in the \code{\link{formula}} by their corresponding functional name
                 }
  \item{data}{an optional \code{\link{data.frame}}, \code{\link{list}} or \code{\link{matrix}} containing the variables of the model, 
	      if not provided the \code{\link{formula}} will be evaluated from the \code{\link{environment}} of the initial \code{\link{call}}}
  \item{weights}{an optional vector of model weights that should be used for estimation}
  \item{family}{specify the distribution for the model, e.g. \code{family="binomial"} to estimate a model of the binomial family with a probit link function}
  \item{iter}{the number of iterations of the gibbs sampler}
  \item{burnin}{the burnin period which will be the cutoff from the single draws of the gibbs sampling}
  \item{thinning}{a parameter which specifies how many samples will be saved during runtime of the gibbs sampler, e.g. if thinning is
                  set to 10, only every 10th draw of the gibbs sampler will be saved. This argument might be important if the data
                  set is very large}
  \item{hyperasigma}{the hyperprior \code{a} of the conjugate prior inverse gamma distribution of the variance \eqn{\sigma^2}{sigma^2}}
  \item{hyperbsigma}{the hyperprior \code{b} of the conjugate prior inverse gamma distribution of the variance \eqn{\sigma^2}{sigma^2}}
  \item{hyperatau}{the hyperprior \code{a} of the conjugate prior inverse gamma distribution of the variance parameter \eqn{\tau^2}{tau^2}}
  \item{hyperbtau}{the hyperprior \code{b} of the conjugate prior inverse gamma distribution of the variance parameter \eqn{\tau^2}{tau^2}}
  \item{trace}{if set to \code{TRUE} every 1000th iteration will be printed in the \R window, however this information will only be
               printed if the buffered output option in the \R window is switched off}
  \item{dots}{print a dot after every 1000th iteration of the sampler}
  \item{speed}{speed up sampling if set to \code{TRUE}}
  \item{\dots}{unused}
}

\details{
  The function estimates a Bayesian hierarchical generalized structured additive model. This is implemented to account for 
  problems with centering of the nonlinear smooth functions. To estimate parameters and functions during the gibbs sampling efficiently, the design matrices
  are transformed in advance, such that the structure of the cross product and penalty matrices will be of diagonal form. This fact reduces estimation time tremendously, 
  since drawing from the full conditionals is not based on multivariate distributions. For details of the transformation functions also see \code{\link{trans.design}}

  Random effects may be estimated within a hierarchical structure of the gibbs sampler, it is possible to model up to three different stages hierarchically, 
  for more details also see \code{\link{r}}
}

\value{
  After the gibbs sampling a list with draws of the coefficients, fitted values, among others will be supplied and stored in the gibbs object.
  The usual generic extractor functions \code{\link{fitted}}, \code{\link{coef}} and \code{\link{residuals}} can be applied to reach estimation results. The user can also
  type \code{names(gibbsobject)} to see all names of the returned objects, the function \code{\link{gibbsObject}} gives a detailed description of the returned
  objects. The information stored in the gibbs objects is of major use of the generic functions \code{\link{summary.gibbs}}, \code{\link{plot.gibbs}} and \code{\link{print.gibbs}}
}

\references{...See the paper...}

\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}}

\seealso{\code{\link{s}}, \code{\link{te}}, \code{\link{sps}}, \code{\link{ste}}, \code{\link{ma}}, \code{\link{mrf}}, \code{\link{m}}, \code{\link{r}}}

\examples{
# estimate a bayesian p-spline
set.seed(121)
n<-500
x1<-runif(n,0,1)

# generate some nonlinear function
f1<-function(x){
		X<-0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
		X<-X-mean(X)}

# generate data 
y<-10+f1(x1)+rnorm(n,sd=1.2)

# estimate model
b<-gibbs(y~s(x1,bs="ps"))
summary(b)
plot(b)

# same with univariate kriging
b<-gibbs(y~s(x1,bs="ma"))
summary(b)
plot(b)


# a 2d example
set.seed(121)
n<-500
f2<-function(x,z){
	sin(x)*cos(z)}

# generate data and plot the true function
x<-runif(n,-3,3)
z<-runif(n,6,10)
xs<-seq(-3,3,length=30)
zs<-seq(6,10,length=30)

pr<-data.frame(x=rep(xs,30),z=rep(zs,rep(30,30)))
truth<-matrix(f2(pr$x,pr$z),30,30)

y<-10+f2(x,z)+rnorm(n,sd=0.2)

# estimate model with tensor product splines
# and with 2d matern splines
# and compare it with gam()
# note multiple penalties is not supported by gibbs()
b1<-gibbs(y~te(x,z,mp=FALSE))
b2<-gibbs(y~s(x,z,bs="ma"))
b3<-gam(y~te(x,z,mp=FALSE))
b4<-gam(y~s(x,z,bs="ma"))
summary(b1)
summary(b2)
summary(b3)
summary(b4)
par(mfrow=c(3,2),mar=c(0,0,2,0))
plot(b1,main="Tensor gibbs() estimate")
plot(b2,main="2d matern gibbs() estimate")
vis.gam(b3,main="Tensor gam() estimate",theta=40,phi=40)
vis.gam(b4,main="2d matern gam() estimate",theta=40,phi=40)
persp(xs,zs,truth,theta=40,phi=40);title("Truth")
par(mfrow=c(1,1))


# estimate simple random effects model
set.seed(111)
n<-1000
     
# index vector for random effects
N<-100
id<-rep(1:N,n/N)
     
# create some iid normal random 
# effects
re<-rnorm(N,sd=0.6)

# add some covariate effects
x1<-sort(runif(n,-3,3))
f1<-sin(x1)

y<-10+f1+re[id]+rnorm(n,sd=0.6)
     
# estimate model
b<-gibbs(y~-1+s(x1)+r(id,~1))
summary(b)

# extract effects and plot 
# against true random effects
# together with true effect
# of x1
par(mfrow=c(1,2))
res<-b$fout[[2]]$effects[,2]
ylim<-c(min(c(res,re)),max(c(res,re)))
plot(res~re,ylim=ylim,xlim=ylim)
lines(c(re)~re,col=2)
plot(b)
lines(f1~x1,col=2)


# estimate model with multilevel effects
set.seed(333)
n<-1000
levels<-6
id<-round(runif(n,1,levels))
x1<-runif(n,-3,3)

M<-rep(0,n)
M[id==1]<-cos(x1[id==1])
M[id==2]<-0.55*x1[id==2]
M[id==3]<-sin(x1[id==3])
M[id==4]<-0.1*(x1[id==4])^2
M[id==5]<-(-0.15)*(x1[id==5])^2
M[id==6]<--0.33*(x1[id==6])

y<-10+M+rnorm(n,sd=1.2)

# estimate model
b<-gibbs(y~m(id,s(x1)))
summary(b)
plot(b,resid=TRUE,cex=0.1)

# same is possible with
id<-as.factor(id)
b<-gibbs(y~s(x1,by=id,xt=list(type="m")))
plot(b)


# estimate a model with regional information
# load a simple sample map
data(polyloc)

# create the corresponding neighborhoodmatrix
neighbor<-neighbormatrix(polyloc)

# generate some values for different regions
n<-1000
x1<-runif(n,0,1)
ind<-round(runif(n,1,length(polyloc)))
m<-nlevels(as.factor(ind))

loc<-rep(0,n)
truth<-sort(runif(m,0,3))
for(i in 1:m)
	loc[ind==i]<-truth[i]
mloc<-mean(loc)
loc<-loc-mloc
truth<-truth-mloc

f2<-function(x){X<-0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
                X<-X-mean(X)}

# generate data
y<-10+f2(x1)+loc+rnorm(n,sd=1.2)	

# estimate model
b<-gibbs(y~s(x1)+s(ind,bs="mrf",xt=list(neighbor)))
summary(b)

par(mfrow=c(1,2))
plot(b,which=1)
plot(b,which=2,map=polyloc,main="Estimated regional effects",
     colored=TRUE,names=TRUE,pal="heat")

# True values
truth


# experiment with weights to account for heteroskedasticity
set.seed(3)
n<-1000
x<-sort(runif(n,3,6))
sdx<-0.2*x
f<-function(x){x<-2*x;x<-x-mean(x)}
e<-rnorm(n,sd=sdx)
y<-f(x)+e

b<-gibbs(y~x)
summary(b)
plot(b,resid=TRUE,cex=0.1)

res<-residuals(b)
ressqrt<-res^2

resmodel<-gibbs(ressqrt~x)
summary(resmodel)
plot(resmodel,resid=TRUE,cex=0.1)

W<-exp(fitted(resmodel))

newmodel<-gibbs(y~x,weights=W)
summary(newmodel)
par(mfrow=c(1,2))
plot(b,resid=TRUE,main="Unweighted model",cex=0.1)
plot(newmodel,resid=TRUE,main="Weighted model",cex=0.1)


# a probit example
n<-1000
x1<-runif(n,-3,3)

z<-sin(x1)+rnorm(n)
y<-rep(0,n)
y[z>0]<-1

b<-gibbs(y~s(x1),family="binomial")
summary(b)
plot(b)


# example with monotonicity 
# constraint on P-spline
set.seed(121)
n<-500
x1<-sort(runif(n,0,1))

fx1<-log(x1)
fx1[(n-2):n]<-fx1[1:2]
y<-10+fx1+rnorm(n,sd=1.2)

# estimate model
b1<-gibbs(y~s(x1,bs="ps"))
summary(b1)
plot(b1,resid=TRUE)

# now impose constraint
b2<-gibbs(y~s(x1,bs="ps",xt=list(constr=1)))
summary(b2)
par(mfrow=c(1,2))
plot(b1,resid=TRUE)
plot(b2,resid=TRUE)
}

\keyword{models}
\keyword{smooth} 
\keyword{regression} 
