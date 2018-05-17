\name{drawmap}
\alias{drawmap}

\title{Drawing Geographical Information}

\description{
Visualises variables that are spatially aligned according to a given map object. Each of the regions in a map will be coloured accoring to the value of the variable.
}
\usage{
drawmap(data, map, regionvar=2, plotvar=3, limits, cols="hcl", nrcolors=100, 
        swapcolors=FALSE, pcat=FALSE,
        hcl.par=list(h=c(120, 0), c=60, l=c(45,60), power=1.2), 
        hsv.par=list(s=1, v=1), legend=TRUE, drawnames=FALSE, cex.names=0.7, 
        cex.legend=0.7, mar.min=2, density=15, ...)}

\arguments{
  \item{data}{Either the name of a file or a data frame containing the variables to be visualised. If missing, the map will be visualised without
              superposition of any further information}
  \item{map}{Map object containing the required boundary information (as obtained by a call to \code{\link{read.bnd}}}
  \item{regionvar}{Defines the variable specifying the geographical regions. Either the name of a variable in \code{data} or the index of the corresponding column.}
  \item{plotvar}{Defines the variable that should be visualised. Either the name of a variable in \code{data} or the index of the corresponding column.}
  \item{limits}{Restricts (or extends) the coloring scheme to a range of values.}
  \item{cols}{Color scheme to be employed. Could be either a vector of colors or one out of the following pre-defined schemes: \code{hcl}, \code{hsv}, \code{grey}}
  \item{nrcolors}{Number of colors (only meaningful when using one of the pre-defined colour schemes).}
  \item{swapcolors}{Reverse the order of colors (works also with user-specified colours but will be most usefule with the pre-defined schemes).}
  \item{pcat}{Option for the visualisation of posterior probabilities. In this case, a three-colour scheme representing significantly positive, insignificant and significantly negative values.}
  \item{hcl.par}{Parameters for the \code{hcl} colour scheme (see the documentation of \code{diverge_hcl} in package \code{vcd} for details).}
  \item{hsv.par}{Parameters for the \code{hsv} colour scheme (see the documentation of \code{\link{hsv}} for details).}
  \item{legend}{Should a legend be added to the figure?}
  \item{drawnames}{Adds the name of each region as a text label to the plot. In most cases the result will be confusing but may be useful when checking the validity of a map.}
  \item{cex.names}{Magnification to be used for the names (if \code{drawnames=TRUE}).}
  \item{cex.legend}{Magnification to be used for the legend.}
  \item{mar.min}{Controls the definition of boundaries. Could be either \code{NULL} for individual settings of \code{mar} or a value which defines \code{mar} as follows: The boundaries will be calculated according to the height to width ratio of the map with minimal boundary \code{mar.min}.}
  \item{density}{Regions without data will be visualised with diagonal stripes. \code{density} defines how dense the stripes should be.}
  \item{...}{Further arguments to be passed to the plot calls that visualise the region boundaries (probably not useful at all).}
}

\author{Felix Heinzl, Thomas Kneib, Andreas Brezger}

\seealso{\code{\link{read.bnd}}}

\examples{
germany <- read.bnd(system.file("examples/germany.bnd", package="BayesX"))
drawmap(map=germany)
drawmap(map=germany, drawnames=TRUE)

res <- read.table(system.file("examples/spatial_f_regions_spatial.res", 
                              package="BayesX"), header=TRUE)

drawmap(res, map=germany)
drawmap(res, map=germany, limits=c(-2,4))
drawmap(res, map=germany, regionvar="regions", plotvar="pmed")
drawmap(res, map=germany, legend=FALSE)
drawmap(res, map=germany, legend=FALSE, main="spatial effect")

drawmap(res, map=germany, cols="hsv")
drawmap(res, map=germany, swapcolors=TRUE, cols="hsv")
drawmap(res, map=germany, cols="grey")
drawmap(res, map=germany,
        cols=c('darkgreen','green','yellow','orange','red','darkred'))

drawmap(res, map=germany, pcat=TRUE, cols="hcl")
drawmap(res, map=germany, pcat=TRUE, cols="hsv")
drawmap(res, map=germany, pcat=TRUE, cols="grey")

drawmap(res, map=germany, nrcolors=10, cols="hcl")
drawmap(res, map=germany, nrcolors=10, cols="hsv")
drawmap(res, map=germany, nrcolors=10, cols="grey")

drawmap(res, map=germany, cols="hcl",
        hcl.par=list(h=c(0,120), c=60, l=c(45,90), power=1.2))
drawmap(res, map=germany, cols="hcl",
        hcl.par=list(h=c(300,120), c=60, l=c(45,90), power=1.2))
drawmap(res, map=germany, cols="hcl",
        hcl.par=list(h=c(40,260), c=60, l=c(45,90), power=1.2))
drawmap(res, map=germany, cols="hsv", hsv.par=list(s=0.7, v=0.7))
}

\keyword{hplot}

