\name{plotsample}
\alias{plotsample}
\alias{plotsample.coda}

\title{Plotting Sampling Paths}

\description{
Plots sampling paths obtained with MCMC-sampling in BayesX
}

\usage{
plotsample(data, ask = TRUE, ...)
plotsample.coda(data, ask = TRUE, ...)
}

\arguments{
  \item{data}{Either the name of a file or a data frame containing the sample.}
  \item{ask}{\code{plotsample} will plot separate sampling paths for each parameter. If \code{ask=TRUE}, 
            the user will be prompted before showing the next plot.}
  \item{\dots}{Further parameters to be passed to the internal call of \code{plot} such as \code{ylim}, etc.}
}

\details{
\code{plotsample} simply plots sampling paths while \code{plotsampe.coda} makes use of the plotting facilities available in package \code{coda}.
}

\author{Felix Heinzl, Andreas Brezger, Thomas Kneib}

\seealso{\code{\link{drawmap}},\code{\link{plotautocor}},\code{\link{plotnonp}},\code{\link{plotsurf}},}

\examples{
res <- read.table(system.file("examples/nonparametric_f_x_pspline_sample.raw",
                              package="BayesX"), header=TRUE)
plotsample(res)
}

\keyword{hplot}
