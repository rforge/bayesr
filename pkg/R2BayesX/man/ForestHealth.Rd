\name{ForestHealth}
\alias{ForestHealth}
\alias{forest}

\title{Forest Health Data}

\description{
  The data set consists of 16 variables and 1796 observations on forest health to identify potential
  factors influencing the health status of trees and therefore the vital status of the forest. In 
  addition to covariates characterizing a tree and its stand, the exact locations of the trees are
  known. The interest is on detecting temporal and spatial trends while accounting for further 
  covariate effects in a flexible manner.
}

\usage{data("ForestHealth")}

\format{
A data frame containing 1796 observations on 16 variables.
\describe{
  \item{age:}{age of stands in years.}
  \item{year:}{age of stands in years.}
  \item{elevation:}{elevation (meters above sea level).}
  \item{inclination:}{slope inclination in percent.}
  \item{soil:}{soil layer depth in cm.}
  \item{ph:}{soil pH at 0-2cm depth.}
  \item{canopy:}{forest canopy density in percent.}
  \item{stand:}{stand type with categories `deciduous' and `mixed'.}
  \item{fertilized:}{fertilization applied with categories `male' and `female'.}
  \item{humus:}{humus layer thickness with categories 1 == `low' to 9 == `high'.}
  \item{moisture:}{soil moisture level with categories 1 == `moderately dry', 2 == `moderately
     moits' and 3 == `moist or temporarily wet'.}
  \item{saturation:}{base saturation with categories 1 == `low' to 4 == `high'.}
}
}

\source{
  \url{http://www.stat.uni-muenchen.de/~bayesx/}.
}

\references{
Kandala, N. B., Lang, S., Klasen, S., Fahrmeir, L. (2001): 
Semiparametric Analysis of the Socio-Demographic and Spatial Determinants of Undernutrition
in Two African Countries. Research in Official Statistics, \bold{1}, 81--100.
}

\seealso{\code{\link{bayesx}}}

\examples{
\dontrun{
## load zambia data and map
data("ZambiaNutrition")
data("ZambiaBnd")

## estimate model
b <- bayesx(stunting ~ memployment + education + urban + gender + 
  s(bmi, bs = "ps") + s(agechild, bs = "ps") +
  s(district, bs = "mrf", xt = list(map = ZambiaBnd)) + r(district),
  iter = 12000, burnin = 2000, step = 10, data = ZambiaNutrition)

summary(b)

## plot smooth effects
plot(b, term = c("s(bmi)", "s(agechild)", "s(district)"), map = ZambiaBnd)

## for more examples 
demo("zambia")
}
}

\keyword{datasets}
