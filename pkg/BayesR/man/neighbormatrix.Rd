\name{neighbormatrix}
\alias{neighbormatrix}

\title{Create a Neighbormatrix from Map}

\description{
The function takes a given map object, which consists of a \code{\link{list}} of \code{\link{polygon}}s, and constructs from the polygon edge and node or centroid information the corresponding neighborhoodmatrix.
}
 
\usage{neighbormatrix(map, type = 1, scale = NULL)}

\arguments{\item{map}{The map object}
           \item{type}{\code{type = 1} will calculate neighbors from sharing borders of the regions, if set to \code{type = 2}, centroids of the \code{\link{polygon}}s 			will be calculcated and used for neighborhood information. 
                        Note that this option is much faster and more suitable for large \code{map} objects. If \code{type = 3} Delaunay Triangulation is used for 				calculating neighbors, also see \code{\link[geometry]{delaunayn}}}
           \item{scale}{if supplied, this must be a number between 0 and 1, this option is only suitable if centroids are calculated, i.e. distances of every
                        centroid are calculated and neighbor centroids are controlled by argument \code{scale}. E.g. \code{scale = 0.1} means that all 
                        centroids which have a distance of less than 10 percent of the maximum distance of two centroids are identified as neighbors}
          }

\value{A neighborhoodmatrix is returned}

\author{Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}}

\seealso{\code{\link{mrf}}}

\examples{
# open map file of nigeria
data(nigeria)

# create the neighbormatrix
nm<-neighbormatrix(nigeria)
nm

# see on map a calculated neighbor dependency
drawmap(nigeria,names=TRUE,main="Neighbors of region '109'")

# draw grey color in region '109'
polygon(nigeria$'109',col=rgb(0.01,0.01,0.01,0.2))

# draw green color into the neighbor polygons of region 109
names<-names(nigeria)
nregions<-length(names)
namesind<-1:length(nregions)
region109<-nm[9,]
for(i in 1:nregions)
	if(region109[i] > 0)
		{
		neighbor<-names[i]
		poly<-eval(parse(text=paste("nigeria$'",neighbor,"'",sep="")))
		polygon(poly,col=rgb(0.01,1,0.01,0.2))
		}
		
# calculate a neighborhoodmatrix with centroids
nm<-neighbormatrix(nigeria,type=2,scale=0.3)

drawmap(nigeria,names=TRUE,main="Neighbors of region '109'")

polygon(nigeria$'109',col=rgb(0.01,0.01,0.01,0.2))

names<-names(nigeria)
nregions<-length(names)
namesind<-1:length(nregions)
region109<-nm[9,]
for(i in 1:nregions)
	if(region109[i] > 0)
		{
		neighbor<-names[i]
		poly<-eval(parse(text=paste("nigeria$'",neighbor,"'",sep="")))
		polygon(poly,col=rgb(0.01,1,0.01,0.2))
		} 
		
nm<-neighbormatrix(nigeria,type=3)
}

\keyword{regression}
\keyword{models}
