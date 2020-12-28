\name{FlashAustria}
\alias{FlashAustria}
\alias{FlashAustriaTrain}
\alias{FlashAustriaEval}
\alias{FlashAustriaCase}

\title{Flash Counts in Austria 2010--2018}

\description{
Lightning discharges counted on 32 x 32 sq-km grid over Austria, Europe,
for the years 2010--2018 and the month May--August, detected by the ground-based
Austrian Lightning Detection and Information System (ALDIS).
A selection of atmospheric parameters from the fifth ECMWF re-analysis (ERA5)
interpolated nearest neighbour method to the grid of the counted lightning
discharges.
The temporal resolution of the data is hourly.
}

\details{
Only a random subset of the data given positive counts
is provided for training and evaluation:
\code{FlashAustriaTrain} and \code{FlashAustriaEval} contain 12000 and 6000
observations, respectively. 
The following prediction cases are provided: 2019-09-15 06 UTC, 2001-09-15 17 UTC, 2001-09-15 23 UTC and
2001-09-16 13 UTC.
}

\format{
Three objects \code{FlashAustriaTrain}, \code{FlashAustriaEval} and
\code{FlashAustriaCase} are added to the environment. \code{FlashAustriaTrain} and
\code{FlashAustriaEval} are \code{data.frame}s. \code{FlashAustriaCase} is an
\code{"sf"} object, containing a geometry column providing the 32 x 32 sq-km grid
as polygons.

\describe{
  \item{counts}{ALDIS counts of cloud-to-ground lightning discharges.}
  \item{hcc}{High cloud cover. The proportion of a grid box covered by cloud occurring in the high levels of the troposhere (approx. 6 km and above).}
  \item{str}{Surface net thermal radiation. Thermal radiation refers to radiation emitted by the atmoshere, clouds and the surface of the Earth. str is the difference between downward and upward thermal radiation at the surface of the Earth. The units are joules per square metre as the radiation is accumulated over a period of one hour.}
  \item{tciw}{Total column cloud ice water. The amount of ice contained within clouds in a column extending from the surface of the Earth to the top of the atmosphere.}
  \item{q_prof_PC1}{This variable has been derived from the vertical profile of specific humidity (labelled q at the ECMWF) by principal component analysis. This is the first principal component.}
  \item{sqrt_cape}{The square root of convective available potential energy. This is an indication of the (in)stability of the atmosphere.}
}
}

\references{
To cite the flash count data (observations), please use:

Schulz, W., Cummins, K., Diendorfer, G., and Dorninger, M. (2005).
Cloud-to-Ground Lightning in Austria: A 10-Year Study Using Data from a Lightning Location System.
J. Geophys. Res., 110, D09101.
\url{https://doi.org/10.1029/2004JD005332}

Schulz, W. and Diendorfer, G. and Pedeboy, S. and Poelman, D. R. (2016).
The European Lightning Location System EUCLID---Part 1: Performance Analysis and Validation.
Nat. Hazards Earth Syst. Sci., 16, 595--605.
\url{https://doi.org/10.5194/nhess-16-595-2016}

To cite the reanalyses data (covariates), please use:

Copernicus Climate Change Service (C3S) (2017).
ERA5: Fifth Generation of ECMWF Atmospheric Reanalyses of the Global Climate.
Copernicus Climate Change Service Climate Data Store (CDS), date of access: June 2019.
\url{https://cds.climate.copernicus.eu/cdsapp#!/home}
}

\examples{
data("FlashAustria", package = "FlashAustria")
head(FlashAustriaTrain)
nrow(FlashAustriaTrain)
hist(FlashAustriaTrain$counts,
  breaks = c(seq(.5, 50.5), Inf),
  xlim = c(0, 50),
  col = "lightgray"
)

plot(counts ~ log(sqrt_cape), FlashAustriaTrain)

nrow(FlashAustriaEval)

if(require("sf")) {
  class(FlashAustriaCase)
  head(FlashAustriaCase)
  plot(st_geometry(FlashAustriaCase))
}

}

\keyword{datasets}

