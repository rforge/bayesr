\name{ste}
\alias{ste}

\title{Create Design and Penalty Matrix for a Tensor Smooth}

\description{
Given two covariates, this function creates the design and penalty matrix for a tensor smooth. The number of knots and the order of the penalty
can be chosen by the user}

\usage{
ste(x, z, by = NULL, degree = c(3, 3), knots = c(20, 20), 
    orderpenalty = c(2, 2), lambda = NULL, center = TRUE, dgts = TRUE,
    map = NULL, ind = NULL)}

\arguments{
  \item{x}{the first covariate the smooth is a function of}
  \item{z}{the second covariate the smooth is a function of}
  \item{by}{this argument specifies if the smooth should be estimated in terms of a varying effect, e.g. the single effect for one factor, then the user
            must supply a dummy variable vector which relates the observations of the covariate to the right factor}
  \item{degree}{the degree of the marginal polynomials}
  \item{knots}{the number of equdistant knots for the marginal basis functions design matrices}
  \item{orderpenalty}{the order of the penalty for the marginal penalty matrices}
  \item{lambda}{if supplied, the value of the hyperparameter which controls for the penalties of the estimated coefficients}
  \item{center}{specifies if the smooth should be centered or not}
  \item{dgts}{feature to reduce the number of total observations of \code{x} and therefore estimation time, 
              if set to \code{TRUE} only different observations in \code{x} will be used for estimation, 
              if set to e.g. \code{dgts = 4} all observations in \code{x} will be rounded to the fourth decimal place and afterwards only different observations in \code{x} will be used for estimation,
              if set to e.g. \code{dgts = c("bin",50)} 50 equidistant design points which cover the range of \code{x} will be calculated and used for etsimation}
  \item{map, ind}{if a \code{map} and the corresponding index is supplied instead of arguments \code{x} and \code{z}, a geospline with centroids as design points of each \code{\link{polygon}} of the \code{map} will be estimated, also see \code{\link{mrf}}}}


 
\details{
The function first constructs the marginal smooth constructors of the two covariates as described in \code{\link{sps}}. In the next step, the tensor
product of the marginal basis design matrices are constructed, which build the new design matrix for the tensor smooth. The new penalty matrix from the 
marginal penalty matrices will be constructed by the sum of kronecker products of the marginal penalty matrices with its identity matrices}

\value{
The following objects are returned by the function
  \item{basis}{the basis function design matrix}
  \item{K}{the penalty matrix}
  \item{x}{the original first covariate}
  \item{z}{the original second covariate}
  \item{degree}{the degree of the marginal polynomials chosen}
  \item{knots}{the number of equidistant knots of the marginals}
  \item{orderpenalty}{the order of the penalties of the marginals}
  \item{names}{the names of the varaibles from the function call}
  \item{lambda}{if supplied, the value of the hyperparameter which will be used during estimation}
  \item{\dots}{ }}

\references{ ~put references to the literature/web site here ~ }

\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}}

\seealso{\code{\link{gibbs}}, \code{\link{ps}}}

\examples{
# estimate a model with tensor splines
set.seed(121)
n<-500

f2<-function(x,z,sx=0.3,sz=0.4)
  { 
  x<-x*20
  (pi**sx*sz)*(1.2*exp(-(x-0.2)^2/sx^2-(z-0.3)^2/sz^2)+
  0.8*exp(-(x-0.7)^2/sx^2-(z-0.8)^2/sz^2))
  }

# generate data and plot the true function
x<-runif(n)/20
z<-runif(n)
xs<-seq(0,1,length=30)/20
zs<-seq(0,1,length=30) 

pr<-data.frame(x=rep(xs,30),z=rep(zs,rep(30,30))) 
truth<-matrix(f2(pr$x,pr$z),30,30) 

y<-10+f2(x,z)+rnorm(n,sd=0.2)

# estimate model
model<-gibbs(y~s(x,z,bs="ste"))
summary(model)

par(mfrow=c(1,2))
plot(model,main="estimate")
persp(xs,zs,truth,theta=40,phi=40)
title("truth")

# plot some basis functions
t<-te(pr$x,pr$z,knots=5)
basis<-t$basis
b<-rep(0,nrow(basis))
check<-5
for(j in 1:ncol(basis))
  if(j==check)
    {
    b<-b+basis[,j]
    check<-check+5
    }
b<-matrix(b,30,30) 
persp(xs,zs,b,theta=40,phi=40,
      xlab="x",ylab="z",zlab="b(x,z)")


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
b<-gibbs(y~s(x1)+s(id,germany,bs="ste",xt=list(geo=TRUE)))
summary(m)

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
