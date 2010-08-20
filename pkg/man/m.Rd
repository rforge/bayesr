\name{m}
\alias{m}

\title{Create Design and Penalty Matrices for Multilevel Effects}

\description{
The function \code{m()} sets up the design and penalty matrices for the \code{\link{gibbs}} sampler. 
Furthermore any method supplied by this package and \code{\link{mgcv }}, i.e. \code{\link{s}}, \code{\link{te}} or simple 
linear effects can be used to construct the corresponding design matrices. 
If for example an effect is supposed to vary within 3 different levels, then this function will supply 3 different design and penalty matrices which 
cover only the domain of the covariate(s) of the certain levels
}
 
\usage{
m(fac, method = NULL, sp = NULL, center = TRUE, hypident = FALSE)
}

\arguments{
  \item{fac}{the factor indicator variable, this variable indicates whether an observation belongs to some certain level or not}
  \item{method}{for linear effects the user can type a formula like argument, e.g. the multilevel linear effect constructor for variable \code{x1} and 
                \code{x2} can be optained by simply typing \code{method = x1 + x2}. If the multilevel effect constructor is of nonlinear type, the
                argument \code{method} has to be specified in terms of the smooth function constructor, i.e. if the multilevel constructor uses
                a \code{\link{ps}} for variable \code{x1}, then \code{method = s(x1)}}
  \item{sp}{the smoothing parameter which controls for the penalties of the estimated coefficients of the smooth terms. If supplied, this argument may be a vector with the
                same length as existing levels. Different hyperparameters can be chosen for the single levels}
  \item{center}{should smooth functions specified in \code{method} be centered or not if set to \code{FALSE}, a vector of specifications may also be supplied as
                in argument \code{lambda}} 
  \item{hypident}{if set to \code{TRUE} only one smoothing parameter will be estimated for all levels}
}

\details{
Supposing there are different levels available, which have different effects each on the response, this function splits up the domain of some given
covariate to the levels specified by the indicator variable fac. E.g. for a 3 level indicator variable the covariate will be split up to the corresponding 3 levels.
In the next step design and penalty matrices, if necessary,  are constructed for every level by the method specified in \code{method} 
}

\value{
Different types of design matrices specified by the argument \code{method} will be returned, as well as a transformator matrix which is used in the 
\code{\link{gibbs}} sampler function, beside other specific objects that will be used during estimation
}

\references{ ~put references to the literature/web site here ~ }

\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}}

\seealso{\code{\link{gibbs}}}

\examples{
# estimate model with multilevel effects
set.seed(111)
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


# now with some 2d effects
set.seed(111)
n<-1000
levels<-3
id<-round(runif(n,1,levels))
                
x<-runif(n,-3,3)
z<-runif(n,6,10)

M<-rep(0,n)
M[id==1]<-sin(x[id==1])*sin(z[id==1])
M[id==2]<-cos(x[id==2])*sin(z[id==2])
M[id==3]<-sin(x[id==3])*0.5*z[id==3]

y<-10+M+rnorm(n,sd=0.6)

b<-gibbs(y~m(id,te(x,z,mp=FALSE)))
summary(b)
par(mar=c(0,0,2,0))
plot(b,resid=TRUE,cex=0.1)

# now only plot level 2
par(mfrow=c(1,1),mar=c(0,0,0,0))
plot(b,which=c(1,2),resid=TRUE,cex=0.1)


# estimate a model with multilevel regional information
n<-3000
data(polyloc)

neighbor<-neighbormatrix(polyloc)

id<-round(runif(n,1,length(polyloc)))
nl<-nlevels(as.factor(ind))

loc<-rep(0,n)

truth1<-c(3,2,2.5,12,13,14,11,15,5,10)
truth1<-truth1 - mean(truth1)
truth2<-c(6,8,2,7.3,1,4,9.2,1.4,3,6.7)
truth2<-truth2 - mean(truth2)
truth3<-c(1.22,9,3.78,4.34,1.3,3,10.03,17,1.23,12.3)
truth3<-truth3 - mean(truth3)
truth<-cbind(truth1,truth2,truth3)

idM<-rep(0,n)
for(i in 1:n)
        {
        idM[i]<-round(runif(1,1,3))
        loc[i]<-truth[id[i],idM[i]]
        }

# generate data
y<-10+loc+rnorm(n,sd=1.2)     

b<-gibbs(y~m(idM,s(id,bs="mrf",xt=list(neighbor))))
summary(b)
par(mar=c(0,0,1,0))
plot(b,map=polyloc,main=c("Level 1","Level 2","Level 3"), 
     colored=FALSE,names=TRUE)

truth1
truth2
truth3

\dontrun{
# example of estimating a number of functions
# for different levels with the same
# smoothing parameter each
set.seed(111)
n<-6000
id<-round(runif(n,1,80))
levels<-nlevels(as.factor(ind))
x1<-runif(n,-3,3)

fx1<-rep(0,n)
err <- rnorm(levels)

# add some variation to the sinus of x1
for(j in 1:levels)
	fx1[id==j] <- err[j]*sin(x1[id==j])

y<-fx1+rnorm(n,sd=1.2)

# estimation may take some time...
b<-gibbs(y~m(id,s(x1),hypident=TRUE))
summary(b)
plot(b,resid=TRUE,cex=0.1)

# same as
id<-as.factor(id)
b<-gibbs(y~s(x1,by=id,xt=list(type="m",hypident=TRUE)))}
}

\keyword{models} 
\keyword{smooth} 
\keyword{regression} 

