\name{plotmap}
\alias{plotmap}

\title{
  Plot Maps
}

\description{
  The function takes a \code{\link[base]{list}} polygons and draws the corresponding map. Different 
  colors for each polygon can be used.
}
 
\usage{
plotmap(map, x = NULL, id = NULL, c.select = NULL, legend = TRUE,
  missing = TRUE, swap = FALSE, range = NULL, names = FALSE,
  values = FALSE, col = NULL, ncol = 100, breaks = NULL,
  cex.legend = 1, cex.names = 1, cex.values = cex.names, digits = 2L,
  mar.min = 2, add = FALSE, interp = FALSE, grid = 200,
  land.only = FALSE, extrap = FALSE, outside = FALSE, type = "akima",
  linear = FALSE, k = 40, p.pch = 15, p.cex = 1, shift = NULL,
  trans = NULL, ...)
}

\arguments{
  \item{map}{The map to be plotted, the map object must be a \code{\link[base]{list}} of matrices 
    with first column indicating the \code{x} coordinate and second column the \code{y} coordinate
    each, also see \code{\link[graphics]{polygon}}.}
  \item{x}{A matrix or data frame with two columns, first column indicates the region and 
    second column the the values which will define the background colors of the polygons, e.g., 
    fitted values from \code{\link{bamlss}}. More columns are possible, e.g. quantiles, which can 
    accessed with argument \code{se}.}
  \item{id}{If argument \code{x} is a vector, argument \code{id} should contain a character vector 
    of the same length of \code{x} with entries indicating the polygon the \eqn{i}-th value of 
    \code{x} belongs to, i.e. \code{id} must contain the same names as polygon names in \code{map}.}
  \item{c.select}{Select the column of the data in \code{x} which should be used for plotting, may
    be an integer or character with the corresponding column name.}
  \item{legend}{If set to \code{TRUE}, a legend will be shown.}
  \item{missing}{should polygons be plotted for which no data is available in \code{x}?}
  \item{swap}{If set to \code{TRUE}, colors will be represented in reverse order.}
  \item{range}{Specify the range of values in \code{x} which should enter the plot, e.g. only 
    values between -2 and 2 are of interest then \code{range = c(-2, 2)}.}
  \item{names}{If set to \code{TRUE} the name for each polygon will also be plotted at the centroids
    of the corresponding polygons.}
  \item{values}{If set to \code{TRUE} the corresponding values for each polygon will also be plotted 
    at the centroids of the polygons.}
  \item{col}{The color of the surface, may also be a function, e.g. 
    \code{col.surface = heat.colors}.}
  \item{ncol}{The number of different colors that should be generated if \code{col} is a function.}
  \item{breaks}{A set of breakpoints for the colors: must give one more breakpoint than 
    \code{ncol}.}
  \item{cex.legend}{Text size of the numbers in the legend.}
  \item{cex.names}{Text size of the names if \code{names = TRUE}.}
  \item{cex.values}{Text size of the names if \code{values = TRUE}.}
  \item{digits}{Specifies the legend decimal places.}
  \item{mar.min}{Controls the definition of boundaries. Could be either \code{NULL} for individual 
    settings of \code{mar} or a value which defines \code{mar} as follows: The boundaries will be 
    calculated according to the height to width ratio of the map with minimal boundary 
    \code{mar.min}.}
  \item{add}{If set to \code{TRUE}, the map will be added to an existing plot.}
  \item{interp}{Logical, should the values provided in argument \code{x} be interpolated to obtain
    a smooth colored map.}
  \item{grid}{Integer, defines the number of grid cells to be used for interpolation.}
  \item{land.only}{If set to \code{TRUE}, only interpoltated pixels that cover land are drawn,
    see also function \code{\link[maps]{map.where}}.}
  \item{extrap}{Logical, should interpolations be computed outside the observation area
    (i.e., extrapolated)?}
  \item{outside}{Logical, should interpolated values outside the boundaries of the map be plotted.}
  \item{type}{Character, which type of interpolation metjod should be used. The default is
    \code{type = "akima"}, see function \code{\link[akima]{interp}}. The two other options are
    \code{type = "mba"}, which calls function \code{\link[MBA]{mba.surf}} of package \pkg{MBA}, or
    \code{type = "mgcv"}, which uses a spatial smoother withing package \pkg{mgcv} for interpolation.
    The last option is definitely the slowest, since a full regression model needs to be estimated.}
  \item{linear}{Logical, should linear interpolation be used withing function
     \code{\link[akima]{interp}}?}
  \item{k}{Integer, the number of basis functions to be used to compute the interpolated surface
    when \code{type = "mgcv"}.}
  \item{p.pch}{Numeric, the point size of the grid cells when using interpolation.}
  \item{p.cex}{Numeric, the size of the grid cell points whein using interpolation.}
  \item{shift}{Numeric constant to be added to the smooth before plotting.}
  \item{trans}{Function to be applied to the smooth before plotting, e.g., to transform the
    plot to the response scale.}
  \item{\dots}{Parameters to be passed to \code{\link{colorlegend}} and others, e.g. change the 
    \code{border} of the polygons and \code{density}, see \code{\link[graphics]{polygon}}. Please 
    see the exmaples.}
}

\note{
  Function \code{plotmap} can use the \pkg{akima} package to construct smooth interpolated
  surfaces, therefore, package \pkg{akima} needs to be installed. The \pkg{akima} package has an ACM
  license that restricts applications to non-commercial usage, see

  \url{http://www.acm.org/publications/policies/softwarecrnotice}
  
  Function \code{plotmap} prints a note refering to the ACM licence. This note can be supressed by
  setting

  \code{options("use.akima" = TRUE)} 
}

\seealso{
  \code{\link{plot2d}}, \code{\link{plot3d}},
  \code{\link{sliceplot}}, \code{\link{plotblock}}.
}

\examples{
## Example from mgcv ?mrf.
## Load Columbus Ohio crime data (see ?columbus for details and credits).
data("columb", package = "mgcv")
data("columb.polys", package = "mgcv")

## Plot the map.
plotmap(columb.polys)

## Plot aggregated data.
a <- with(columb, aggregate(crime,
  by = list("district" = district), FUN = mean))

plotmap(columb.polys, x = a$x, id = a$district)
plotmap(columb.polys, x = a$x, id = a$district,
  pos = "topleft")
plotmap(columb.polys, x = a$x, id = a$district,
  pos = "topleft", side.legend = 2)
plotmap(columb.polys, x = a$x, id = a$district,
  pos = "topleft", side.legend = 2, side.ticks = 2)
plotmap(columb.polys, x = a$x, id = a$district,
  pos = "topleft", side.legend = 2, side.ticks = 2,
  col = heat_hcl, swap = TRUE)
plotmap(columb.polys, x = a$x, id = a$district,
  pos = "topleft", side.legend = 2, side.ticks = 2,
  col = heat_hcl, swap = TRUE, range = c(10, 50))
plotmap(columb.polys, x = a$x, id = a$district,
  pos = "topleft", side.legend = 2, side.ticks = 2,
  col = heat_hcl(5), swap = TRUE, range = c(10, 50),
  lrange = c(0, 60))

\dontrun{## Estimate model.
xt <- list("polys" = columb.polys)
b <- bamlss(crime ~ s(district,bs="mrf",xt=xt) + s(income), data = columb)

## Plot effects.
plot(b, map = columb.polys, pos = "topleft", digits = 0)
}}

\keyword{hplot}

