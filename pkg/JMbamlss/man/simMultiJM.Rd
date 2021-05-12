\name{simMultiJM}
\alias{simMultiJM}

\title{
  %%  ~~function to do ... ~~
}

\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
simMultiJM(nsub = 300, times = seq(0, 120, 1), probmiss = 0.75,
  nmark = 2, M = 6, ncovar = 2,
  lambda = function(t, x) { 1.4 * log((t + 10)/1000) - 1.5 },
  gamma = function(x) { 0.3 * x[, 1] },
  alpha = rep(list(function(t, x) { 0.3 + 0 * t}), nmark),
  mu = rep(list(function(t, x) { 1.25 + 0.6 * sin(x[, 2]) + (-0.01) * t }), nmark),
  sigma = function(t, x) { 0.3 + 0 * t + I(x$marker == "m2") * 0.2 },
  tmax = NULL, seed = NULL, mfpc_args = list(type = "split",
  eFunType = "Poly", ignoreDeg = NULL, eValType = "linear", eValScale = 1),
  full = FALSE, file = NULL)
}

\arguments{
  \item{nsub}{
%%     ~~Describe \code{nsub} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{probmiss}{
%%     ~~Describe \code{probmiss} here~~
}
  \item{nmark}{
%%     ~~Describe \code{nmark} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{ncovar}{
%%     ~~Describe \code{ncovar} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{tmax}{
%%     ~~Describe \code{tmax} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{mfpc_args}{
%%     ~~Describe \code{mfpc_args} here~~
}
  \item{full}{
%%     ~~Describe \code{full} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\references{
%% ~put references to the literature/web site here ~
}

\author{
%%  ~~who you are~~
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  ## Add example here!
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

