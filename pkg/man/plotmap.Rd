\name{plotmap}
\alias{plotmap}

\title{Plot Maps}

\description{
  The function takes a \code{\link{list}} of \code{\link{polygon}}s and draws the corresponding map. Different colors for each polygon can be used.
  The function is mainly used internal by the plotting function \code{\link{plot.gibbs}} for a \code{\link{gibbs}} object
}
 
\usage{
plotmap(map, names = FALSE, values = NULL, colored = FALSE, range = NULL, 
        pal = "rainbow", legend = TRUE, scale = 0.2, nrc = 100, xlab = NULL, ylab = NULL, 
	zlab = NULL, main = NULL, dgts = 4, cex = 1, lpos = c(0.2,0), p3d = FALSE, theta = 0, 
	phi = 25, pcat = NULL, border = NULL, density = NULL, resid = FALSE, ...)
}

\arguments{
  \item{map}{the map to be plotted, the map object must be a \code{\link{list}} of \code{\link{polygon}}s}
  \item{names}{if set to \code{TRUE} the regionames will also be plotted}
  \item{values}{the values, e.g. fitted values, which will define different colors for the polygons}
  \item{colored}{if set to \code{TRUE} a colored map will be drawn}
  \item{range}{specifiy a range of \code{values} which should enter the plot, e.g. only \code{values} between 3 and 5 are of interest
               then \code{range = c(3,5)}}
  \item{pal}{specifies which palette should be used for the colors, options are \code{"rainbow","heat","terrain","topo","cm","hcl"}}
  \item{legend}{if set to \code{TRUE} a legend will be shown}
  \item{scale}{scales the size of the map, negative scales will increase the view}
  \item{nrc}{the number of different colors of the chosen palette}
  \item{xlab}{a label for the x axis}
  \item{ylab}{a label for the y axis}
  \item{zlab}{a label for the z axis if \code{p3d} is set to \code{TRUE}}
  \item{main}{a main title for the plot}
  \item{dgts}{rounds the legend values to the specified number of decimal places}
  \item{cex}{if names is set to \code{TRUE} the size of the text can be specified, also see \code{\link{plot.default}}}
  \item{lpos}{specifies the position of the legend, useful values range between 0 and 1, first element specifies the position on the x-axis and second on the y-axis}
  \item{p3d}{if argument \code{map} is supplied, this option will plot the map in 3d using the functionality of \code{\link{persp}}}
  \item{theta, phi}{angles defining the viewing direction. theta gives the azimuthal direction and phi the colatitude}
  \item{border}{choose the borders of the polygons, see \code{\link{polygon}}}
  \item{density}{also see \code{\link{polygon}}}
  \item{resid}{if set to \code{TRUE}, partial resiuals will be plotted if available}
  \item{\dots}{other graphical parameters to be passed on, e.g. specifify the parameters for the \code{hcl} colour scheme by providing a list named like 
	       \code{hcl.par=list(h=c(130,25), c=100, l=c(90,70))} (see the documentation of \code{diverge_hcl} in package \code{vcd} for details).
		It is also possible to reverse the order of colors by setting \code{swap = T}}
}

\value{
  A map from the \code{\link{polygon}} \code{\link{list}} in \code{map} will be drawn
}

\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}}

\seealso{\code{\link{plot.gibbs}}, \code{\link{readbndfile}}}

\examples{
# load a sample map
data(polyloc)

# plot the map 
plotmap(polyloc,main="Example of a plain map")

# add some values to the corresponding polygon areas
# note that the first column in matrix val contains
# the region identification index
val<-cbind(as.integer(names(polyloc)),runif(length(polyloc),0,100))
plotmap(polyloc,values=val, 
        main="Example of a map with black and white shaded regions")
        
# now only plot values for some certain regions
samps<-sample(val[,1],4)
nval<-val[samps,]
plotmap(polyloc,values=nval)

# play with scale and legend position
plotmap(polyloc,values=nval,scale=0,cex=1,dgts=2,lpos=c(0.6,0),density=15)

# choose a color palette and print region names
plotmap(polyloc,values=val, 
        main="Example of a map with color palette heat", 
        colored=TRUE,pal="heat",names=TRUE)

# load and plot a map from germany
data(germany)
plotmap(germany,main="Map of Germany")

# draw the map in 3d
plotmap(germany,main="Map of Germany",p3d=TRUE)

# add some colors
val<-cbind(as.integer(names(germany)),sort(runif(length(germany),0,100)))
plotmap(germany,values=val,colored=TRUE)

# remove polygon borders
plotmap(germany,values=val,colored=TRUE,border=NA)

# use hcl colors
plotmap(germany,values=val,colored=TRUE,border=NA,pal="hcl")

# swap color
plotmap(germany,values=val,colored=TRUE,border=NA,pal="hcl",swap=TRUE)
}

\keyword{aplot} 
