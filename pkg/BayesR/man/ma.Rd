\name{ma}
\alias{ma}

\title{Create Design and Penalty Matrices for Unavariate or Bivariate Matern Kriging Methods}
\description{
 This function generates the design and penalty matrices which are required for the \code{\link{gibbs}} routine within function \code{\link{s}} 
 and \code{\link{smooth.construct}}. The correlation functions in the kriging method are from the matern family correlation functions.  
 An algorithm which chooses optimal knot locations is used for univariate and bivariate kriging. 
 This algorithm is adopted from the package \code{\link[fields]{fields}}, see also \code{\link[fields]{cover.design}}}
 
\usage{
ma(x, z, by = NULL, knots = 20, c = 11.75, kappa = 2, 
   loc = NULL, lambda = NULL, center = TRUE, dgts = TRUE,
   map = NULL, ind = NULL)

# and within gibbs()
s(..., bs="ma", xt=NULL}

\arguments{
  \item{x}{the covariate that the matern smooth is a function of}
  \item{z}{second possible covariate that the matern smooth is a function of}
  \item{by}{this argument specifies if the smooth should be estimated in terms of a varying effect, e.g. the single effect for one factor, then the user
            must supply a dummy variable vector which relates the observations of the covariate to the right factor}
  \item{knots}{the number of knots for the smooth}
  \item{c}{parameter which controls the range of one basis function}
  \item{kappa}{parameter which controls for the type of the correlation function, e.g. if \code{kappa = 0.5}, then an exponential correlation function is applied}
  \item{loc}{if existing, a certain number of knot locations can be chosen by the user which will be used for estimation}
  \item{lambda}{the hyperparameter which controls for the penalties on the estimated coefficients, if supplied only this parameter will be used}
  \item{center}{specifies if the smooth should be centered or not}
  \item{dgts}{feature to reduce the number of total observations of \code{x} and therefore estimation time, 
              if set to \code{TRUE} only different observations in \code{x} will be used for estimation, 
              if set to e.g. \code{dgts = 4} all observations in \code{x} will be rounded to the fourth decimal place and afterwards only different observations 
	      in \code{x} will be used for estimation, if set to e.g. \code{dgts = c("bin",50)} 50 equidistant design points which cover the range of \code{x} 
	      will be calculated and used for estimation}
  \item{map, ind}{if a \code{map} and the corresponding index is supplied instead of arguments \code{x} and \code{z}, a geospline with centroids as design 
		  points of each \code{\link{polygon}} of the \code{map} will be estimated, also see \code{\link{mrf}}}
  \item{xt}{here, specifify all arguments to be modified from default specifications within ma(), i.e. set the \code{loc} argeument with 
	    \code{xt=list(loc=myloc))}}}

\value{
  Function ma() returns the design and penalty matrix to set up a \code{\link{gibbs}} model, the most important objects returned are
  \item{basis}{the design matrix}
  \item{by}{if supplied, a variable the smooth is varying by}
  \item{K}{the corresponding penalty matrix}
  \item{dim}{identifier variable, e.g. \code{dim=1} indicates univariate kriging}
  \item{x}{the original covariate}
  \item{z}{if provided, the second original covariate}
  \item{kappa}{the value chosen for parameter kappa}
  \item{names}{the names of the covariate(s)}
  \item{lambda}{if supplied, the value of the hyperparameter}
  \item{knots}{the nmuber of knots chosen in \code{knots}}
  \item{loc}{the knot locations}
  \item{\dots}{ }}


\references{ ~put references to the literature/web site here ~ }

\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}}

\seealso{\code{\link{gibbs}}}
\examples{
# create a matern basis funtion design matrix 
n<-100
x<-sort(runif(n,0,10))
b<-ma(x,knots=NULL)
summary(b)

# plotting of basis functions
basis<-b$basis
matplot(x=x,y=basis,type="l",pch=1,ylab="b(x)")

# little variation
b<-ma(x,kappa=0.5,knots=NULL)
basis<-b$basis
matplot(x=x,y=basis,type="l",pch=1,ylab="b(x)")

# now with two covariates and 40 knots
# using knot point searching with cover.design()
n<-900
x<-runif(n,-3,3)
z<-runif(n,-3,3)
basis<-ma(x,z,knots=40) 

# plot knot locations choosen from cover.design()
# in red color togther with all design points
plot(z~x)
points(basis$loc,col="red",pch=16)

# now plot the corresponding basis functions 
# in a little movie with design points and
# knot locations
grid<-40
xs<-seq(min(x),max(x),length=grid)
zs<-seq(min(z),max(z),length=grid)
pr<-data.frame(x=rep(xs,grid),z=rep(zs,rep(grid,grid))) 
b2d<-ma(pr$x,pr$z,loc=basis$loc)
for(j in 1:ncol(basis$basis))
  {
  bb<-matrix(b2d$basis[,j],grid,grid) 

  res<-persp(x=xs,y=zs,z=bb,theta=40,phi=40,
             ticktype="detailed",xlab="x",ylab="z",
             zlab="")
  text(-0.5,0, parse(text=paste("B[",j,"](x,z)",sep="")),srt=-73)
  points(trans3d(x=x,y=z,z=0,pmat=res),col=1)
  points(trans3d(basis$loc[,1],basis$loc[,2],0,pmat=res),col=3,pch=10)
  points(trans3d(basis$loc[j,1],basis$loc[j,2],0,pmat=res),col=2,pch=16)
  Sys.sleep(0.3)
  }


# a 1d example
set.seed(121)
n<-500
x1<-runif(n,0,1)

# generate some nonlinear function
f1<-function(x){
		X<-0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
		X<-X-mean(X)}

# generate data 
y<-10+f1(x1)+rnorm(n,sd=1.2)

# estimate model with univariate kriging
b<-gibbs(y~s(x1,bs="ma"))
summary(b)
plot(b)

# change parameters for c and kappa 
# but use same knot locations
# and change c with m[1] <- 1
# and kappa with m[2] <- 1
loc <- attr(fitted(m,type=2)[[1]],"smooth.specs")$xt$loc
b<-gibbs(y~s(x1,bs="ma",m=c(1,1),xt=list(loc=loc)))
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

b<-gibbs(y~s(x,z,bs="ma"))
summary(b)

par(mfrow=c(1,2))
plot(b,main="2d matern gibbs() estimate")
persp(xs,zs,truth,theta=40,phi=40);title("Truth")


\dontrun{
# now a geospline example
set.seed(111)

# first simulate some geographical data
data(germany)
N<-length(germany)
names(germany)<-1:N
n<-N*5
id<-rep(1:N,n/N)
f2<-sort(runif(N,0,3))
f2<-f2-mean(f2) # center for better comparison
f2<-f2[id]

# generate some univariate data
x1<-runif(n,-3,3)
f1<-sin(x1)

y<-10+f1+f2+rnorm(n,sd=0.6)

# estimate model with geospline kriging
b<-gibbs(y~s(x1)+s(id,germany,bs="ma",k=(N-1),xt=list(geo=TRUE)))
summary(b)

f<-fitted(b,type=2)[[2]]
cp<-geospline.centroid.data
fs<-cbind(cp[order(cp[,1]),3],f[!duplicated(f[,1]),3])

par(mfrow=c(1,3))
plot(germany,values=fs,colored=TRUE,pal="hcl",main="Estimate",dgts=3)
plot(germany,values=cbind(1:N,f2),colored=TRUE,pal="hcl",main="Truth",dgts=3)
plot(f[!duplicated(f[,1]),3]~f2[order(cp[id,1])][!duplicated(f[,1])],
     xlab="True effects",ylab="Fitted")
lines(c(f2)~f2,col=2)}
}

\keyword{models} 
\keyword{smooth} 
\keyword{regression} 
