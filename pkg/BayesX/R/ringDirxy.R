## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* campus *.* lmu *.*
## de] Project: BayesX Time-stamp: [ringDirxy.R] by DSB Sam 28/02/2009 19:18 (GMT)
## on daniel@puc-home> Description: Private R code replacing
## maptools:::.ringDirxy.  History: 28/02/2009 file creation
.ringDirxy <- function(xy) {
  xy <- as.matrix(xy)
  nPoints <- nrow(xy)
  
  ## checks
  stopifnot(is.numeric(xy), identical(ncol(xy), as.integer(2)), nPoints >= 3)
  
  ## index vector from the 2nd to the 2nd but last point
  inds <- seq_len(nPoints - 2) + 1
  
  ## the start point and the middle and tail matrix
  start <- xy[1, ]  # here the drop of dims is OK
  middle <- xy[inds, , drop = FALSE]  # here prevent drop of dims if there are only 3 points
  tail <- xy[inds + 1, , drop = FALSE]
  
  ## compute twice the signed areas
  areas <- (middle[, 1] - start[1]) * (tail[, 2] - start[2]) - (tail[, 1] - start[1]) * 
    (middle[, 2] - start[2])
  ## and sum up to total twice signed area of the polygon
  total <- sum(areas)
  
  ## the sign then gives the direction
  if (total > 0) {
    return(as.integer(-1))  # counter-clockwise
  } else {
    return(as.integer(1))  # clockwise
  }
}

