\name{plot.gibbs}
\alias{plot.gibbs}
\alias{plot.sm.gibbs}
\alias{plot.mrf.gibbs}
\alias{plot.random.gibbs}
\alias{plot.linear.gibbs}

\title{Default Gibbs Plotting}

\description{
The generic functions take fitted \code{\link{gibbs}} objects and plots a user selected term of the model, also parametric model components can be plotted
}
 
\usage{
\method{plot}{gibbs}(x, which = 1, resid = FALSE, jit = TRUE, const = FALSE, diagnostics = FALSE, 
                     acf = FALSE, xlab = NULL, ylab = NULL, zlab = NULL, main = NULL, 
                     colored = TRUE, col = FALSE, lwdc = 1, lwdconf = 0, grid = 30, theta = 40, 
                     phi = 40, image = FALSE, map = NULL, names = FALSE, values = NULL, range = NULL, 
                     pal = "heat", legend = TRUE, scale = 0.2, nrc = 100, dgts = 2, lpos = c(0.2,0),
                     cex = 1, ask = FALSE, byplots = FALSE, p3d = FALSE, pcat, border = NULL, ...) 
     
\method{plot}{linear.gibbs}(x, resid = FALSE, jit = TRUE, const = FALSE, 
                            diagnostics = FALSE, acf = FALSE, xlab = NULL, ylab = NULL, 
                            main = NULL, colored = TRUE, col = FALSE, lwdc = 1, 
                            lwdconf = 0, range = NULL, cex = 1, ask = FALSE, ...) 
			                   
\method{plot}{mrf.gibbs}(x, resid = FALSE, map = NULL, names = FALSE, values = NULL, 
                         colored = TRUE, col = FALSE, lwdc = 1, lwdconf = 0, range = NULL, 
                         pal = "heat", legend = TRUE, scale = 0, nrc = 100, xlab = NULL, 
                         ylab = NULL, zlab = NULL, main = NULL, dgts = 4, cex = 1, lpos = c(0,0), 
                         const = FALSE, diagnostics = FALSE, acf = FALSE, p3d = FALSE, theta = 0, 
                         phi = 25, pcat, border = NULL, ...)                         
			                   
\method{plot}{random.gibbs}(x, resid = FALSE, const = FALSE, diagnostics = FALSE, acf = FALSE, 
                            xlab = NULL, ylab = NULL, zlab = NULL, main = NULL, colored = TRUE, 
                            col = FALSE, lwdc = 1, lwdconf = 0, map = NULL, names = FALSE, 
                            range = NULL, pal = "heat", legend = TRUE, scale = 0.2, nrc = 100, 
                            dgts = 2, lpos = c(0.2,0), cex = 1, ask = FALSE, p3d = FALSE, theta = 0, 
                            phi = 25, pcat, border = NULL, ...)

\method{plot}{sm.gibbs}(x, resid = FALSE, jit = TRUE, const = FALSE, diagnostics = FALSE, 
                        acf = FALSE, xlab = NULL, ylab = NULL, zlab = NULL, main = NULL, 
			colored = TRUE, col = FALSE, lwdc = 1, lwdconf = 0, grid = 30, 
			theta = 40, phi = 40, image = FALSE, cex = 1, ask = FALSE, ...)
}

\arguments{
  \item{x}{a fitted \code{\link{gibbs}} object}
  \item{which}{which term should be plotted, e.g. the model formula looks like \code{y ~ x + ps(z)}, if the user wants to plot the smooth function
               estimated by the pspline, then \code{which = 2} must be selected, for \code{diagnostics} plots, \code{which} has the same meaning}
  \item{resid}{if set to \code{TRUE} partial residuals are plotted}
  \item{jit}{if set to \code{TRUE}, a \code{\link{rug}} will be added to the plot, only possible if \code{resid = FALSE}}
  \item{const}{if set to \code{TRUE}, the effect with all other terms set to the \code{mean} will be plotted}
  \item{diagnostics}{if diagnostic plots should be supplied, if set to 1, residual diagnostic plots will be displayed, if set to 2, 
                     a sampling path of the variance, if the term selected in \code{which} is a smooth, is provided, if set to 3 sampling paths of the coefficicients 
                     of the corresponding term will be plotted, if set to 4 the sampling path of the coefficient of the intercept is drawn}
  \item{acf}{use \code{acf = TRUE} if \code{diagnostics = 2}, which will plot autocorrelation functions of the sampled coefficients, also see \code{\link{acf}}}
  \item{xlab}{a label for the \code{x} axis}
  \item{ylab}{a label for the \code{y} axis}
  \item{zlab}{a label for the \code{z} axis for 2D plots}
  \item{main}{a main title}
  \item{colored}{if set to \code{TRUE}, the space between the confidence bands will be filled with color}
  \item{col}{specifies the color that should be printed between the confidence bands, e.g. a red area can be achieved by 
                 \code{col = "red"}, etc.}
  \item{lwdc}{specify the linewidth of the main effect}
  \item{lwdconf}{specify the linewidth of the confidence bands}
  \item{grid}{for 2d plotting, specify the number of data points in the 2 dimensional plot}
  \item{theta, phi}{angles defining the viewing direction. theta gives the azimuthal direction and phi the colatitude}
  \item{image}{if set to \code{TRUE}, a color image (see \code{\link{image}}) instead of 2 dimensional plot will be supplied together with contours}
  \item{map}{if the \code{\link{gibbs}} object term is of type \code{\link{mrf}} then a map list object can be supplied to plot the map with the
             estimated regional effects, also see \code{\link{drawmap}}. This is also possible for \code{\link{ra}}ndom effects plots}
  \item{names}{if a map object is supplied, the regionnames are plotted if \code{names} is set to \code{TRUE}, also see \code{\link{drawmap}}}
  \item{values}{the values, e.g. fitted values, which will define different colors for the polygons}
  \item{range}{specifies the range of a map plot, also see \code{\link{drawmap}}. If the term to be plotted is a \code{\link{factor}}, argument \code{range} may specify the range of the boxes drawn for the single factor effect}
  \item{pal}{chose a color palette for a map plot, see \code{link{drawmap}}}
  \item{legend}{if set to \code{TRUE} a legend will be supplied in a map plot, see \code{\link{drawmap}}}
  \item{scale}{scales the map plot, see \code{\link{drawmap}}}
  \item{nrc}{the number of different colors in a map plot, see \code{\link{drawmap}}}
  \item{dgts}{if a map is supplied, rounds the legend values to the specified number of decimal places, see \code{\link{drawmap}}}
  \item{lpos}{specifies the position of the legend, useful values range between 0 and 1, first element specifies the position on the x-axis and second on the y-axis}
  \item{cex}{specifiy the text size of regionnames in a map plot, also see \code{\link{drawmap}}, and for other usage see \code{\link{plot.default}}}
  \item{ask}{should multiple plots be supplied one after another, see also \code{\link{par}}}
  \item{byplots}{if argument \code{by} in function \code{\link{ra}} is provided, this argument specifies if plots for \code{by} should be drawn}
  \item{p3d}{if argument \code{map} is supplied, this option will plot the map in 3d using the functionality of \code{\link{persp}}}
  \item{pcat}{if spatial estimates should be plotted within a map object, \code{pcat = 95} will plot the 95 percent posterior propabilies only, \code{pcat = 80} 80 percent posterior probailites. Note that a similar result is achieved by setting \code{nrc = 3} and e.g. \code{range = c(-1,1)}}
  \item{border}{choose the borders of the polygons, see \code{\link{polygon}}}
  \item{\dots}{other graphical parameters to be passed on}
}

\details{
This function automatically identifies the type of the term the user selects by the \code{which} argument. If for an example a 2 dimensional plot is required,
a plot from the function \code{\link{persp}} or \code{\link{image}} will be called. It os also possible to extract a fitted term object and plot it, this is
also possible with \code{\link{predict.gibbs}}
}

\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}}

\seealso{\code{\link{gibbs}}, \code{\link{predict.gibbs}}}

\examples{
# generate some covariates
set.seed(111)
n<-500
x<-runif(n,-3,3)
z<-runif(n,-3,3)
w<-runif(n,0,6)
fac<-round(runif(n,1,10))
Z<-diag(nlevels(fac))[fac,]

# generate some nonlinear function 1d and 2d
f1<-function(x){sin(x)}
f2<-function(x,z){cos(x)*sin(z)}

# some values for the factor covariate
betas<-c(2.67,5,6,3,4,2,6,7,9,7.5)[fac]
fac<-as.factor(fac)

# generate data 
y<-10+f1(x)+f2(z,w)+betas+rnorm(n,sd=0.6)

# estimate model
m<-gibbs(y~s(x)+te(z,w,mp=FALSE)+fac)

# plot pspline term
plot(m,which=1)

# with residuals
plot(m,which=1,resid=TRUE)

# plot tensor term
plot(m,which=2)

# plot tensor term with residuals
plot(m,which=2,resid=TRUE)

# plot image and contour
plot(m,which=2,image=TRUE)
plot(m,which=2,image=TRUE,contour=TRUE)

# plot factor term
plot(m,3)

# plot factor term with residuals
plot(m,3,resid=TRUE,cex=0.1)

# plot residual dignostics
plot(m,diagnostics=1, 
     main=c("Histogram of residuals",
            "Fitted values vs. residuals",
            "Residual Diagnostics"))

# plot variance sampling path of term s(x)
plot(m,diagnostics=2,which=1)

# plot coefficients sampling paths of term s(x)
plot(m,diagnostics=3,which=1)

# plot the sampling path of the intercept
plot(m,diagnostics=4)

# plot the correlation function of the sampled intercept
plot(m,diagnostics=4,acf=TRUE)

# extract a fitted term object
fx <- fitted(m,type=2)[[1]]

# and plot it, same as plot(m,which=1)
plot(fx)
}

\keyword{aplot}
