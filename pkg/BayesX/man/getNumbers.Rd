\name{getNumbers}
\alias{getNumbers}

\title{Extract numbers from (BayesX log file) strings}
\description{
  This is an internal helper function for extractSamples which extracts
  numbers from (BayesX log file) strings.
}
\usage{
getNumbers(beforeStringsList, stringVector)
}

\arguments{
  \item{beforeStringsList}{The list with the strings standing before the
  numbers which shall be extracted.}
  \item{stringVector}{The vector of strings to be searched for the
  numbers.} 
}
\value{
  Returns a list with the extracted numbers.
}

\author{Daniel Sabanes Bove}

\examples{
## create a nice example
sampleStrings <- c("Second: 385",
                   "  First:  70000 ",                   
                   "asdfkf T: 24      ")

## test the function
extractedNumbers <- BayesX:::getNumbers(beforeStringsList=
                                        list(first="First:",
                                             second="Second:",
                                             third="T:"),
                                        stringVector=sampleStrings)

## and assert that the result is as expected
stopifnot(identical(extractedNumbers,
                    list(first=70000, second=385, third=24)))

## specific test
stopifnot(identical(BayesX:::getNumbers(beforeStringsList=
                                        list(Iterations = "Number of iterations:"),
                                        stringVector=
                                        "  Number of iterations: 70000 "),
                    list(Iterations=70000)))
}

\keyword{internal}

