\name{plotautocor}
\alias{plotautocor}

\title{Computing and Plotting Autocorrelation Functions}
\description{
Computes and plot autocorrelation functions for samples obtained with MCMC in BayesX
}

\usage{
plotautocor(data, ask = TRUE, lag.max=100, ...)
}

\arguments{
  \item{data}{Either the name of a file or a data frame containing the sample.}
  \item{ask}{\code{plotautocor} will plot separate autocorrelation functions for each parameter. If \code{ask=TRUE}, 
            the user will be prompted before showing the next plot.}
  \item{lag.max}{Maximum number of lags to be considered.}
  \item{\dots}{Further parameters to be passed to the internal call of \code{plot} such as \code{ylim}, etc.}
}

\author{Felix Heinzl, Thomas Kneib}

\examples{
res <- read.table(system.file("examples/nonparametric_f_x_pspline_sample.raw", 
                              package="BayesX"), header=TRUE)
plotautocor(res)
plotautocor(res, lag.max=50)
}

\keyword{hplot}
