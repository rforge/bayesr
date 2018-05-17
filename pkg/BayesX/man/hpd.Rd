\name{hpd}
\alias{hpd}
\alias{hpd.coda}

\title{Computing Highest Posterior Density (HPD) Intervals}

\description{
Compute approximate HPD intervals out of MCMC-samples in BayesX
}

\usage{
hpd(data, alpha = 0.05, ...)
hpd.coda(data, alpha = 0.05)
}

\arguments{
  \item{data}{Either the name of a file or a data frame containing the sample.}
  \item{alpha}{A numeric scalar in the interval (0,1) such that 1 - alpha is the
    target probability content of the intervals.. The default is alpha = 0.05.}
  \item{\dots}{Further parameters to be passed to the internal call of \code{optim}
    and \code{integrate}.}
}

\details{
\code{hpd} computes the HPD interval based on a kernel density estimate of the samples.
\code{hpd.coda} computes the HPD interval with the function \code{HPDinterval} available in
  package \code{coda}.
}

\author{Nadja Klein}

\examples{
res <- read.table(system.file("examples/nonparametric_f_x_pspline_sample.raw",
  package="BayesX"), header = TRUE)
hpd(res)
hpd.coda(res)
}

\keyword{MCMC}
