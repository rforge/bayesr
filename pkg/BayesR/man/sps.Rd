\name{ps}
\alias{ps}

\title{Create P-spline Design and Penalty Matrix}

\description{
The function generates the design and penalty matrix for a p-spline smooth from a given covariate. The order of the penalty and the number of equidistant 
knots are chosen by the user
}
 
\usage{
ps(x, by = NULL, degree = 3, knots = 20, orderpenalty = 2, 
   lambda = NULL, center = TRUE, dgts = TRUE)}

\arguments{
  \item{x}{the covariate that the p-spline is a function of}
  \item{by}{this argument specifies if the smooth should be estimated in terms of a varying effect, e.g. the single effect for one factor, then the user
            must supply a dummy variable vector which relates the observations of the covariate to the right factor}
  \item{degree}{the degree of the polynomial}
  \item{knots}{the number of equidistant knots}
  \item{orderpenalty}{the order of the penalty}
  \item{lambda}{the hyperparameter which controls for the penalties on the estimated coefficients, if supplied only this parameter will be used}
  \item{center}{specifies if the smooth should be centered or not}
  \item{dgts}{feature to reduce the number of total observations of \code{x} and therefore estimation time, 
              if set to \code{TRUE} only different observations in \code{x} will be used for estimation, 
              if set to e.g. \code{dgts = 4} all observations in \code{x} will be rounded to the fourth decimal place and afterwards only different observations in \code{x} will be used for estimation,
              if set to e.g. \code{dgts = c("bin",50)} 50 equidistant design points which cover the range of \code{x} will be calculated and used for etsimation}}

\value{
The following arguments are returned by the function
  \item{basis}{the basis function design matrix}
  \item{K}{the penalty matrix}
  \item{x}{the original covariate}
  \item{by}{the dummy variable which specifies the varying effect}
  \item{names}{the name of the covariate and dummy variable if supplied}
  \item{lambda}{if supplied, the value of the hyperparameter}
  \item{degree}{the degree of the pspline}
  \item{knots}{the number of equidistant knots}
  \item{orderpenalty}{the order of the penalty}
  \item{\dots}{ }}
  
\references{ ~put references to the literature/web site here ~ }

\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}}

\seealso{\code{\link{gibbs}}, \code{\link{penalty}}, \code{\link{te}}}

\examples{
# generate data
n<-500
x<-sort(runif(n,0,50))

# generate pspline smooth constructor
mypspline<-ps(x)
summary(mypspline)

# plotting of basis functions
basis<-mypspline$basis
matplot(x=x,y=basis,type="l",pch=1,
        xlab="x",ylab="b(x)")
}

\keyword{models} 
\keyword{smooth} 
\keyword{regression} 
