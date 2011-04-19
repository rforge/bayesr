\name{BayesR-package}
\alias{BayesR-package}
\alias{BayesR}
\alias{bayesr}
\docType{package}
 
\title{Bayesian Estimation of Hierarchical Generalized Structured Additive Regression Models
}

\description{Routine for the estimation of Bayesian Hierarchical Generalized Structured Additive 
             Regression Models with MCMC Gibbs sampling, including optimal smoothing parameter estimation.
}

\details{
\tabular{ll}{
Package: \tab BayesR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-07-22\cr
License: \tab Nik\cr
}
The main estimation function is called \code{\link{gibbs}}, which initializes a \code{\link{gibbs}} sampling algorithm for estimation of models with structured additive predictor with MCMC. 
Like other packages summary statistics and plots are provided for \code{\link{gibbs}} objects (see \code{\link{summary.gibbs}} and \code{\link{plot.gibbs}}). 
Nonlinear smooth functions may be modeled by \code{\link{s}} or \code{\link{te}} terms using the model \code{\link{formula}}e specifications as supplied by package \code{\link{mgcv}}, 
also see \code{\link{smooth.terms}}, where additional bases for e.g. modeling spatial information with markov random fields \code{\link{mrf}}, univariate and bivariate 
\code{\link{ma}}atern Kriging methods, as well as another tensor P-spline constructor for modeling geosplines is provided, see 
\code{\link{ste}}. Additionally, random effects may be modeled in terms of covariates within a hierarchical structure, see \code{\link{r}},
and multilevel effects with the \code{\link{m}} constructor. 

Estimation is thereby based on a transformation of the parameters, which leads to diagonal structures of the design and precision matrices and therefore the sampler only requires to draw random numbers from univariate distributions, only, see \code{\link{trans.design}}. Remarkable speed increase is also achieved by taking into account that the number of different observations of some covariate is usually much smaller than the total number of observations. There are two versions of the sampling algorithm, one is fully implemented in \code{R} code, and for speed advantages there is also an implementation in \code{C} provided.
}

\author{
Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at} and Stefan Lang \email{stefan.lang@uibk.ac.at}

Maintainer: Nikolaus Umlauf \email{nikolaus.umlauf@uibk.ac.at}
}

\references{
Lang, S., Umlauf, N., Kneib, T. and Wechselberger, P. (2010): Multilevel Generalized Structured Additive Regression.
}

\keyword{ package }

\examples{
# see help sites for functions gibbs, s, te, sps, ste, ma, mrf, m, r
}

\keyword{models} 
\keyword{smooth} 
\keyword{regression} 
