\name{ForestHealth}
\alias{ForestHealth}
\alias{forest}

\title{Forest Health Data}

\description{
  The data set consists of 16 variables with 1796 observations on forest health to identify 
  potential factors influencing the health status of trees and therefore the vital status of the 
  forest. In addition to covariates characterizing a tree and its stand, the exact locations of the 
  trees are known. The interest is on detecting temporal and spatial trends while accounting for 
  further covariate effects in a flexible manner.
}

\usage{data("ForestHealth")}

\format{
A data frame containing 1793 observations on 16 variables.
\describe{
  \item{id:}{tree location identification number.}
  \item{year:}{year of census.}
  \item{defoliation:}{percentage of tree defoliation in three ordinal categories, 
    0 == `defoliation < 12.5\%', 1 == `12.5\% <= defoliation < 50\%' and  
    3 == `defoliation >= 50\%'}
  \item{x:}{x-coordinate of the tree location.}
  \item{y:}{y-coordinate of the tree location.}
  \item{age:}{age of stands in years.}
  \item{canopy:}{forest canopy density in percent.}
  \item{inclination:}{slope inclination in percent.}
  \item{elevation:}{elevation (meters above sea level).}
  \item{soil:}{soil layer depth in cm.}
  \item{ph:}{soil pH at 0-2cm depth.}
  \item{moisture:}{soil moisture level with categories 1 == `moderately dry', 2 == `moderately
     moist' and 3 == `moist or temporarily wet'.}
  \item{alkali:}{proportion of base alkali-ions with categories 1 == `low' to 4 == `high'.}
  \item{humus:}{humus layer thickness in cm, categorical coded.}
  \item{stand:}{stand type with categories `deciduous' and `mixed'.}
  \item{fertilized:}{fertilization applied with categories `yes' and `no'.}
}
}

\source{
  \url{http://www.stat.uni-muenchen.de/~kneib/regressionsbuch/}.
}

\references{
Kneib, T. & Fahrmeir, L. (2010): A Space-Time Study on Forest Health.
In: Chandler, R. E. \& Scott, M. (eds.): Statistical Methods for Trend
Detection and Analysis in the Environmental Sciences, Wiley.
}

\seealso{\code{\link{bayesx}}}

\examples{
\dontrun{
## load zambia data and map
data("ForestHealth")
data("BeechBnd")

fm <- bayesx(defoliation ~  stand + fertilized + 
  humus + moisture + alkali + ph + soil + 
  s(age, bs = "ps") + s(inclination, bs = "ps") + 
  s(canopy, bs = "ps") + s(year, bs = "ps") + 
  s(elevation, bs = "ps"),
  family = "cumlogit", method = "REML", data = ForestHealth)

summary(fm)
plot(fm, term = c("s(age)", "s(inclination)", 
  "s(canopy)", "s(year)", "s(elevation)"))
}
}

\keyword{datasets}
