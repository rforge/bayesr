\name{r}
\alias{r}

\title{Create Random Effects Model Matrices}

\description{
This function creates the design and penalty matrices for a random effects model specified in \code{\link{gibbs}}. Given any covariate(s) this
function may use any of the functions described in \code{\link{s}}, \code{\link{te}}, \code{\link{ma}}, \code{\link{mrf}} to set up the terms for hierarchical estimation.
It is possible to set up multiple internal hierarchical structures within the constructor function.
}

\usage{
r(fac, method = NULL, by = NULL)
}
 
\arguments{
  \item{fac}{the factor variable which relates the observations of the random effects to the corresponding response observation}
  \item{method}{if not supplied, only random effects which are not modeled in terms of covariates will be estimated. For linear effects the user can
                type a \code{\link{formula}} type argument, e.g. the random effect constructor for variable \code{x1} and \code{x2} can be obtained by  
                typing \code{method =~ x1 + x2}. If the random effect constructor uses functions of nonlinear type, the argument \code{method} has to be specified
                in terms of the smooth function constructor, i.e. if the \code{r}andom effects constructor uses a \code{\link{ps}} for variable \code{x1},
                then \code{method =~ s(x1) + x2}. It is also possible to model a third stage of random effects, this is initialized by simply adding another
                \code{r} term in the \code{method} \code{\link{formula}}, e.g. \code{method =~ s(x1) + x2 + r(fac2)}, where also additional covariates
                may enter the model within the third stage, see the examples}
  \item{by}{specify a variable for which the slope estimator is modeled as random}
}

\details{
...Write more details...
}

\value{
The following arguments are returned 
  \item{n}{the number of observations from the random effects covariate}
  \item{Z}{a 0/1 incidence design matrix which maps the random effects observation to the repsonse}
  \item{call}{the initial function call}
  \item{fac}{the factor variable which relates the observations of the random effects to the corresponding response observation}
  \item{\dots}{other objects depending on the type of function calls}
}

\references{ ~put references to the literature/web site here ~ }

\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}}

\seealso{\code{\link{gibbs}}}

\examples{
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
m<-gibbs(y~-1+s(x1)+r(id,~1))
summary(m)

# extract effects and plot 
# against true random effects
# together with true effect
# of x1
par(mfrow=c(1,2))
res<-m$fout[[2]]$effects[,2]
ylim<-c(min(c(res,re)),max(c(res,re)))
plot(res~re,ylim=ylim,xlim=ylim)
lines(c(re)~re,col=2)
plot(m)
lines(f1~x1,col=2)

# now estimate a random effects
# model with a stage 2 covariate
# effect
set.seed(333)
n<-1000

# 1st stage
x1<-runif(n,0,6)
f1<-cos(x1)

# 2nd stage
N<-100
x2<-sort(runif(N,-3,3))
id<-rep(1:N,n/N)
f2<-sin(x2)
re<-rnorm(N,sd=0.6)
fre<-f2+re

y<-10+f1+fre[id]+rnorm(n,sd=0.6)

# estimate model
# hierarchical
m1<-gibbs(y~-1+s(x1)+r(id,~1+s(x2)))
summary(m1)

# non hierarchical
X2<-x2[id]
m2<-gibbs(y~-1+s(x1)+s(X2)+r(id,~1))
summary(m2)

# plot stage 2 effects
par(mfrow=c(1,2))
plot(m1,which=c(2,2),resid=TRUE,cex=0.1,main="hierarchical")
lines(f2~x2,col=2)
plot(m2,which=2,resid=TRUE,cex=0.1,main="non hierarchical")
lines(f2~x2,col=2)

# setup a model with 2 random effect stages
set.seed(333)
n<-4000
N<-400
x1<-runif(n,-3,3)
id1<-as.factor(round(runif(n,1,N)))
nr<-nlevels(id1)
x2<-runif(nr,-3,3)

NN<-100
id2<-as.factor(round(runif(nr,1,NN)))
nr2<-nlevels(id2)
x3<-runif(nr2,-2,4)

# some functions for the terms
fr2<-cos(x3)
fr2<-fr2-mean(fr2)
re2<-rnorm(nr2,sd=0.6)
re2<-re2-mean(re2)

fr1<-sin(x2)
fr1<-fr1-mean(fr2)
re1<-rnorm(nr,sd=0.6)
re1<-re1-mean(re1)
random<-c(fr1+re1+c(fr2+re2)[id2])[id1]
y<-cos(x1)+random+rnorm(n,sd=0.6)

# estimate 3 stage hierarchical model
m<-gibbs(y~-1+s(x1)+r(id1,~-1+s(x2)+r(id2,~1+s(x3))),iter=20000,burnin=10000,thinning=10)
summary(m)
par(mfrow=c(1,3))
plot(m,which=1)
lines(cos(x1)[order(x1)]~x1[order(x1)],col=2)
plot(m,which=c(2,1))
lines(fr1[order(x2)]~x2[order(x2)],col=2)
plot(m,which=c(2,2,2))
lines(fr2[order(x3)]~x3[order(x3)],col=2)

# example of estimating a number of functions
# for different levels with the same
# smoothing parameter each
set.seed(111)
n<-6000
ind<-round(runif(n,1,80))
levels<-nlevels(as.factor(ind))
x1<-runif(n,-3,3)

fx1<-rep(0,n)
err <- rnorm(levels)

# add some variation to the sinus of x1
for(j in 1:levels)
	fx1[ind==j] <- err[j]*sin(x1[ind==j])

y<-10+fx1+rnorm(n,sd=1.2)

# not run: estimation will take some time...
# m<-gibbs(y~r(ind,by=~s(x1)))
# summary(m)
# plot(m,resid=TRUE,cex=0.1,byplots=TRUE)

# setup a model with 3 random effect stages
# note, hierarchical estimation only possible if speed = TRUE in gibbs()
set.seed(333)
n<-4000
N<-400
x1<-runif(n,-3,3)
id1<-as.factor(round(runif(n,1,N)))
nr<-nlevels(id1)
r1<-runif(nr,-3,3)

NN<-200
id2<-as.factor(round(runif(nr,1,NN)))
nr2<-nlevels(id2)
r2<-runif(nr2,-2,4)

NNN<-100
id3<-as.factor(round(runif(nr2,1,NNN)))
nr3<-nlevels(id3)
r3<-runif(nr3,6,12)

# some functions for the terms
Z1<-diag(nlevels(id1))[id1,]
Z2<-diag(nlevels(id2))[id2,]
Z3<-diag(nlevels(id3))[id3,]

fr3<-sin(r3)
fr3<-fr3-mean(fr3)
re3<-rnorm(nr3,sd=0.6)
re3<-re3-mean(re3)

fr2<-cos(r2)
fr2<-fr2-mean(fr2)
re2<-rnorm(nr2,sd=0.6)
re2<-re2-mean(re2)

fr1<-sin(r1)
fr1<-fr1-mean(fr2)
re1<-rnorm(nr,sd=0.6)
re1<-re1-mean(re1)
random<-Z1\%*\%(fr1+re1+Z2\%*\%(fr2+re2+Z3\%*\%(fr3+re3)))

y<-10+cos(x1)+random+rnorm(n,sd=0.6)

m1<-gibbs(y~-1+s(x1)+r(id1,~-1+s(r1)+r(id2,~-1+s(r2)+r(id3,~1+s(r3)))))
summary(m1)
par(mfrow=c(2,2))
plot(m1,1)
plot(m1,c(2,1))
plot(m1,c(2,2,1))
plot(m1,c(2,2,2,1))

# now the same none hierarchical
r3<- c(Z1\%*\%(Z2\%*\%(Z3\%*\%r3)))
id3<- c(Z1\%*\%(Z2\%*\%as.integer(id3)))
r2<- c(Z1\%*\%(Z2\%*\%r2))
id2<- c(Z1\%*\%as.integer(id2))
r1<-c(Z1\%*\%r1)

m2<-gibbs(y~-1+s(x1)+s(r1)+s(r2)+s(r3)+r(id,~-1)+r(id2,~-1)+r(id3,~1))
summary(m2)
par(mfrow=c(2,2))
plot(m2,1)
plot(m2,2)
plot(m2,3)
plot(m2,4)
}

\keyword{models} 
\keyword{smooth} 
\keyword{regression} 

