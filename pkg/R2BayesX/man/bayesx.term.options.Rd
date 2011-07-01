\name{bayesx.term.options}
\alias{bayesx.term.options}
\alias{bxopts}

\title{Show BayesX Term Options}

\description{
  BayesX model terms specified using functions \code{\link[mgcv]{s}} and \code{\link{r}} may have 
  additional optional control arguments within argument \code{xt}. The argument passed to \code{xt}
  must be a named list with the specifications of the control parameters. Therefore function 
  \code{bayesx.term.options} displays the possible additional controlling parameters for a 
  particular model term.
}

\usage{
bayesx.term.options(bs = "ps")
}

\arguments{
  \item{bs}{the term specification for which controlling parameters should be shown. Note,
    for random effects specified with function \code{\link{r}}, \code{bs = "random"} will display
    random effects optional arguments.}
}

\details{
  At the moment the following model terms are implemented, for which additional controlling 
  parameters may be specified:
  \itemize{
    \item \code{"lasso"}, \code{"nigmix"}, \code{"ridge"}: shrinkage of fixed effects, defines a 
      shrinkage-prior for the corresponding parameters \eqn{\gamma_j}, \eqn{j = 1, \ldots, q}, 
      \eqn{q \geq 1} of the linear effects \eqn{x_1, \ldots, x_q}. There are three priors possible: 
      ridge-, lasso- and Normal Mixture of inverse Gamma (NMIG)-prior. 
    \item \code{"bl"}: nonlinear baseline effect in hazard regression or multi-state models, defines 
      a P-spline with second order random walk penalty for the parameters of the spline for the 
      log-baseline effect \eqn{log(\lambda(time))}.
    \item \code{"kr"}: kriging with stationary Gaussian random fields. 
    \item \code{"gk"}: geokriging with stationary Gaussian random fields: estimates a stationary 
      Gaussian random field based on the centroids of a map object provided in boundary format, see
      function \code{\link[BayesX]{read.bnd}}, within the \code{xt} argument of function 
      \code{\link[mgcv]{s}}, e.g. \code{xt = list(map = MapBnd)}. 
    \item \code{"gs"}: geosplines based on two-dimensional P-splines with first order random walk 
      penalty, defines a two-dimensional P-spline for the spatial covariate region with a 
      two-dimensional first order random walk penalty for the parameters of the spline. Estimation 
      is based on the coordinates of the centroids of the regions provided by a map object in in 
      boundary format, see function \code{\link[BayesX]{read.bnd}}, within the \code{xt} argument of 
      function \code{\link[mgcv]{s}}, e.g. \code{xt = list(map = MapBnd)}. 
    \item \code{"mrf"}: Markov random fields: defines a Markov random field prior for a spatial 
      covariate, where geographical informations are provided by a map object in boundary or graph 
      file format, see function \code{\link[BayesX]{read.bnd}} and \code{\link[BayesX]{read.gra}}, 
      within the \code{xt} argument of function \code{\link[mgcv]{s}}, e.g. 
      \code{xt = list(map = MapBndorGra)}. 
    \item \code{"ps"}: P-spline with first or second order difference penalty. 
    \item \code{"rw1"}, \code{"rw2"}: zero degree P-splines, defines a zero degree P-spline with 
      first or second order difference penalty. A zero degree P-spline typically estimates for every 
      distinct covariate value in the data set a separate parameter. Usually there is no reason to 
      prefer zero degree P-splines over higher order P-splines. An exception are ordinal covariates 
      or continuous covariates with only a small number of different values. For ordinal covariates 
      higher order P-splines are not meaningful while zero degree P-splines might be an alternative 
      to modeling nonlinear relationships via a dummy approach with completely unrestricted 
      regression parameters.
    \item \code{"season"}: seasonal effect of a time scale. 
    \item \code{"te"}: defines a two-dimensional P-spline based on the tensor product of 
      one-dimensional P-splines with a two-dimensional first order random walk penalty for the 
      parameters of the spline.
  }
}

\author{
  Thomas Kneib, Stefan Lang, Nikolaus Umlauf, Achim Zeileis.
}

\examples{
## show arguments for P-splines
bayesx.term.options(bs = "ps")

## Markov random fields
bayesx.term.options(bs = "mrf")
}

\keyword{regression}
