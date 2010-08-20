\name{predict.gibbs}
\alias{predict.gibbs}

\title{Prediction from fitted Gibbs Model}

\description{
Function to calculate predictions of a fitted \code{\link{gibbs}} model.
}
 
\usage{
\method{predict}{gibbs}(object, which = 1, newdata = NULL, ...)
}

\arguments{
  \item{object}{a \code{\link{gibbs}} object}
  \item{which}{which term should be used, e.g. the model \code{\link{formula}} looks like \code{y ~ x + ps(z)}, if the user wants to calculate the prediction of the smooth function
               estimated by the pspline, then \code{which = 2} must be selected. If the model contains random effects terms modeled hierarchical, the \code{which} argument may be a vector of length depending on the
               depth of the levels specified, i.e. the model contains a third stage, then \code{which} may be specified as \code{which = c(2,2)}, in the case the second term of the model is 
               a random effect and the second argument in the random effects constructor is modeled as another random effect level, e.g. see \code{y ~ x + ra(id1, ~ps(z) + ra(id2))}. If the third stage contains an additional covariate, then
               \code{which} may be specified as \code{which = c(2,2,1)}, e.g. see the \code{\link{formula}} \code{y ~ x + ra(id1, ~ps(z) + ra(id2, ~ps(q)))}. It is also possible to get the predictions of the effects only in this case, simply choose 
               \code{which = c(2,2)}, this is also valid for the second stage, \code{which = 2}. If \code{which = NULL} the predictor is returned}
  \item{newdata}{Supply new data for the term specified in \code{which}, for bivariate functions \code{newdata} must be a \code{\link{matrix}} with 2 columns}
  \item{\dots}{further arguments passed to or from other methods}
}

\value{
This function returns mean and quantiles of the model term predictions estimated with \code{\link{gibbs} given original values of covariates or new data}
}

\references{\dots}

\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}}

\seealso{\code{\link{gibbs}}}

\examples{
# estimate a bayesian p-spline
set.seed(12345)
n<-500
x1<-runif(n,0,1)

# generate some nonlinear function
f1<-function(x){X<-0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
                X<-X-mean(X)}

# generate data 
y<-10+f1(x1)+rnorm(n,sd=1.2)

# estimate model
model<-gibbs(y~ps(x1))
summary(model)

# get predictions
pred<-predict(model,which=1)
plot(pred)

# create some new values for x1
newx1 <- seq(0.1,0.8,length=30)
pred<-predict(model,which=1,newdata=newx1)
plot(pred)

pred
}

\keyword{models} 
\keyword{smooth} 
\keyword{regression} 

