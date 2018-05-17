\name{plotnonp}
\alias{plotnonp}

\title{Plotting Nonparametric Function Estimates}

\description{
Plots nonparametric function estimates obtained from BayesX
}

\usage{
plotnonp(data, x = 2, y = c(3, 4, 5, 7, 8), ylim = NULL, 
         lty = c(1, 2, 3, 2, 3), cols = rep(1, length(y)), month, year, step=12, 
         xlab, ylab, ...)
}

\arguments{
  \item{data}{Either the name of a file or a data frame containing the estimation results.}
  \item{x}{Defines the x-axis in the plot. Either the name of a variable in \code{data} or the index of the corresponding column.}
  \item{y}{Defines the variables to be plotted against \code{x}. May be either a vector of names of variables
           in \code{data} or the corresponding indices. The default choice corresponds to the point estimate plus two
           confidence bands.}
  \item{ylim}{Since \code{plotnonp} plots multiple y-variables, it automatically determines the appropriate
  \code{ylim} to make all curves visible. Argument \code{ylim} allows to override this default behaviour with fixed values.}
  \item{lty}{Vector of line types used for plotting (must have the same length as \code{y}). The default corresponds
        to solid lines for the point estimate and dashed and dotted lines for the confidence bands.}
  \item{cols}{Vector of colors used for plotting (must have the same length as \code{y}). Default are black lines.}
  \item{month, year, step}{Provide specific annotation for plotting estimation results for temporal variables. \code{month}
       and \code{year} define the minimum time point whereas step specifies the type of temporal data with
       \code{step=4}, \code{step=2} and \code{step=1} corresponding to quartely, half yearly and yearly data.}
  \item{xlab, ylab}{\code{plotnonp} constructs default labels that can be overwritten by these arguments}
  \item{\dots}{Further arguments to be passed to the interval call of \code{plot} such as \code{type}, etc.}
}

\author{Felix Heinzl, Andreas Brezger and Thomas Kneib}

\seealso{\code{\link{drawmap}},\code{\link{plotautocor}},\code{\link{plotsample}},\code{\link{plotsurf}}}

\examples{

res <- read.table(system.file("examples/nonparametric_f_x_pspline.res", 
                              package="BayesX"), header=TRUE)
plotnonp(res)
plotnonp(res, x="x")
plotnonp(res, x="x", y="pmean")
plotnonp(res, x="x", y="pmed")
plotnonp(res, x="x", y="pmed", ylim=c(-2,2))
plotnonp(res, x="x", y=c("pmean", "pqu10", "pqu90"), lty=c(1,1,1), 
         col=c("red","blue","blue"))
plotnonp(res, xlab="some variable", ylab="f(some variable)", 
         main="Nonlinear effect of some variable", sub="penalised spline")

res <- read.table(system.file("examples/nonparametric2_f_time_pspline.res", 
                              package="BayesX"), header=TRUE)
plotnonp(res)
plotnonp(res, month=1, year=1980, step=12)

res <- res[1:18,]                                           
plotnonp(res, month=1, year=1980, step=12)
}

\keyword{hplot}
