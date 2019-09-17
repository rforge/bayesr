\name{Golf}
\alias{Golf}

\title{Prices of Used Cars Data}

\description{
  This dataset is taken from the Regression Book and is about prices of used VW Golf cars.
}

\usage{data("Golf")}

\format{
The \code{Golf} data contains the following variables:
\describe{
  \item{price:}{Numeric, sale price in 1000 Euro.}
  \item{age:}{Numeric, Age of the car in month.}
  \item{kilometer:}{Numeric, Kilometer reading in 1000 kilometer.}
  \item{TIA:}{Numeric, Month until the next TIA appointment.}
  \item{abs:}{Factor, Does the car have abs?}
  \item{sunroof:}{Factor, Does the car have a sunroof?}
}
}

\references{
  Fahrmeir, L., Kneib, T., Lang, S. and Marx, B. (2013).
  Regression - Models, Methods and Applications, Springer.
  \url{https://www.uni-goettingen.de/de/books/550514.html}.
}

\examples{
data("Golf")
plot(price ~ age, data = Golf)
}

\keyword{datasets}

