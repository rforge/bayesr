read.gra <- function(file, sorted = FALSE, sep = " ")
{    
  ## read information from the graph file
  gd <- strsplit(readLines(file), sep)

  ## get region names and neighbors
  if(length(gd[[2]]) < 2) {
    ids <- rep(1:3, length = length(gd) - 1)
    rn <- unlist(gd[c(2:length(gd))[ids == 1]])
    ng <- gd[c(2:length(gd))[ids == 3]]
  } else {
    rn <- sapply(gd[-1], function(x) x[1])
    ng <- sapply(gd[-1], function(x) x[-1])
  }
  ng <- lapply(ng, as.integer)

  ## extract the number of districts
  n <- length(rn)
  cat("Note: map consists of", S, "regions\n")

  cat("Creating adjacency matrix ...\n")

  ## create the adjacency matrix
  pmat <- matrix(0, n, n)
  
  ## number of neighbors for each region
  diag(pmat) <- dp <- sapply(ng, length)

  ## loop over the districts
  for(i in 1:n) {
    ## off-diagonal entries in pmat
    ## Note: the indices are zero-based!
    if(length(ng[[i]]))
      pmat[i, ng[[i]]] <- -1
  }

  colnames(pmat) <- rownames(pmat) <- rn

  cat("finished\n")
    
  if(sorted) {
    ## sort districts and adjacency matrix
    if(sum(is.na(as.numeric(districts))) == 0){
      rn <- as.numeric(rn)
      cat("Note: regions sorted by number\n") 
    } else  cat("Note: regions sorted by name\n")
    ord <- order(rn)
    rn <- sort(rn)
        
    pmat.sort <- matrix(0, S, S)

    for(i in 1:n){
      ordi <- ord[i]
      for(j in 1:S){
        ordj <- ord[j]
        pmat.sort[i, j] <- pmat[ordi, ordj]
      } 
    }
        
    pmat <- pmat.sort
  }

  rownames(pmat) <- colnames(pmat) <- rn

  class(pmat) <- "gra"
  return(pmat)
}

