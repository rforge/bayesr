\name{n}
\alias{n}
\alias{smooth.construct.nn.smooth.spec}
\alias{n.weights}

\title{Neural Networks for BAMLSS}

\description{
  This smooth constructor implements single hidden layer neural networks.
}

\usage{
## The neural network smooth constructor.
n(..., k = 10)

## Initialize weights.
n.weights(nodes, k, r = NULL, s = NULL,
  type = c("sigmoid", "gauss", "softplus", "cos", "sin"),
  x = NULL, ...)

## For mgcv.
\method{smooth.construct}{nn.smooth.spec}(object, data, knots, ...)
}

\arguments{
  \item{\dots}{For function \code{n()} a formula of the type \code{~x1+x2+x3} that specifies
    the covariates that should be modeled by the neural network.}
  \item{k}{For function \code{n()}, the number of hidden nodes of the network. Note that one can set
    an argument \code{split = TRUE} to split up the neural network into, e.g., \code{nsplit = 5}
    parts with \code{k} nodes each. For function \code{n.weights()}, argument \code{k}
    is the number of input variables of the network (number of covariates).}
  \item{nodes}{Number of nodes for each layer, i.e., can also be a vector.}
  \item{r, s}{Parameters controlling the shape of the activation functions.}
  \item{type}{The type of activation function that should be used.}
  \item{x}{A scaled covariate matrix, the data will be used to identify the range of the weights.}
  \item{object, data, knots}{See \code{\link[mgcv]{smooth.construct}}.}
}

\value{
  Function \code{n()}, similar to function \code{\link[mgcv]{s}} a simple smooth specification
  object.
}


\seealso{
\code{\link{bamlss}}, \code{\link{predict.bamlss}}, \code{\link{bfit}}, \code{\link{boost}}
}

\examples{
## ... coming soon ...!
}

\keyword{regression}

