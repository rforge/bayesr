\name{simJM}
\alias{simJM}
\alias{rJM}


\title{
  Simulate longitudinal and survival data within the joint model framework
}

\description{
  Simulates longitudinal data with normal error and (Cox-type) survival times
  using the inversion method. The function \code{simJM()} is a wrapper specifying
  all predictors. The wrapper calls  \code{rJM()} which is modified from \code{rSurvtime()}
  of the \R package \code{CoxFlexBoost}. 
}

\usage{
simJM(nsub = 300, times = seq(0, 120, 1), probmiss = 0.75,
      long_setting = "functional", alpha_setting = "nonlinear", 
      dalpha_setting="zero", sigma=0.3, long_df=6, tmax=NULL,    
      seed=NULL, full=TRUE, file = NULL, ...){
}

\arguments{
  \item{nsub}{number of individuals for which longitudinal data and survival times should be simulated.}
  \item{times}{vector of time points at which longitudinal measurements are "sampled".}
  \item{probmiss}{proportion of longitudinal measurements to be set to missing. Used to induce sparsity in the longitudinal measurements.}
  \item{long_setting}{Specification of the longitudinal trajectories of the sampled subjects. Preset specifications are \code{"simple"}, \code{"linear"}, \code{"nonlinear"} and \code{"functional"}. See Details.}
  \item{alpha_setting}{specification of the association between survival and longitudinal. Preset specifications are \code{"simple"}, \code{"linear"}, \code{"nonlinear"} and \code{"nonlinear2"}. See Details.}
  \item{dalpha_setting}{specification of the association between survival and the derivative of the longitudinal. Work in progress.}
  \item{sigma}{standard deviation of the normal error around the true longitudinal measurements.}
  \item{long_df}{number of basis functions from which functional random intercepts are sampled.}
  \item{tmax}{longest possible survival time, observations are censored after that timepoint. 
  Defaults to \code{max(times)}. Should not be specified longer than \code{max(times)} for functional random intercepts.}
  \item{seed}{numeric scalar setting the random seed.}
  \item{full}{logical indicating if only the longitudinal data set should be returned (\code{FALSE}) or also the data for the survival part evaluated on a regular time grid (\code{TRUE}). }
  \item{file}{character. name of the data file the generated data set should be stored into
    (e.g., "simdata.RData") or NULL if the dataset should directly be returned in \R. }
  \item{\dots}{further arguments to be passed to the predictor functions. }

}

\details{
  The function simulates longitudinal data basing on the given specification at given \code{times}. 
  The full hazard is built from all joint model predictors and survival times are sampled using the 
  inversion method (cf. Bender, Augustin, & Blettner, 2005). Additional censoring and missingness 
  can be introduced. The longitudinal information is censored according to the survival information. 
  Pre-specified functions of \code{long_setting} and \code{alpha_setting} are: TODO.
}

\value{
  A \code{data.frame} in long-format consisting of
  \item{id}{subject identifier}
  \item{survtime}{sampled survival time}
  \item{event}{non-censoring indicator}
  \item{cumhaz}{cumulative hazard at the survival time}
  \item{obstime}{time points of longitudinally "observed" measurements}
  \item{y}{longitudinal response}
Further the covariates \code{x1} and \code{x2}, the random effects \code{r1} and \code{r2}
as well as the true predictors \code{alpha}, \code{lambda}, \code{gamma}, \code{mu} and
under construction \code{dmu} and \code{dalpha}.
are returned
  the non-censoring indicator (\code{event}) and further
  covariates \code{x} is returned. If \code{file} is specified, the
  data.frame is additionally stored on the disc.
}

\references{
  Benjamin Hofner (2016). CoxFlexBoost: Boosting Flexible Cox Models
  (with Time-Varying Effects). R package version 0.7-0.

  Ralph Bender and Thomas Augustin and Maria Blettner (2005),
  Generating Survival Times to Simulate Cox Proportional Hazards
  Models.
  \emph{Statistics in Medicine}, \bold{24}, 1713--1723.
}

\seealso{
  \code{\link{jm.bamlss}}, \code{\link{jm.mode}}, \code{\link{jm.mcmc}}, \code{\link{bamlss}}
}

\examples{
## todo
}

\keyword{survival}
\keyword{distribution}

