bnd2gra <- function(map, npoints = 2)
{
    ## check class of argument
    if (! inherits(map, "bnd"))
        stop("argument 'map' is not an object of class 'bnd'")

    ## extract (unique) regions from the polygons
    regions <- unique(names(map))
    nRegions <- length(regions)

    ## initialize return matrix
    retMatrix <- matrix(data=0,
                        nrow=nRegions, ncol=nRegions,
                        dimnames=list(regions, regions))

    ## helper function:
    pointsMatrixToPointsCharVector <- function(pointsMatrix)
    {
        paste(pointsMatrix[, 1L],
              pointsMatrix[, 2L],
              sep="&")
    }

    ## check for neighboring polygons:
    ## process only upper triangular part of the matrix, because it is symmetric.

    ## counter for already processed pairs
    nPairsProcessed <- 0L

    ## how many pairs need to be processed?
    nPairs <- (nRegions * (nRegions - 1L)) / 2L

    ## and at which iterations should a progress message be printed?
    nProgressIterations <- 10L
    progressIterations <- seq(from=0L,
                              to=nPairs,
                              length=nProgressIterations + 1L)[- 1L]

    cat("Start neighbor search ...\n")
    
    pointsMatrix <- NULL

    ## i is the row region number
    for (i in seq_len(nRegions - 1L))
    {
        ## which polygons belong to region number i?
        polyIndsRegion.i <- which(names(map) == regions[i])

        ## make a vector of all points (x, y) in the format "x&y" which belong to this region
        if(length(polyIndsRegion.i)>1)
          {
          pointsMatrix <- map[polyIndsRegion.i[1]]
          for(k in 2:length(polyIndsRegion.i))
            pointsMatrix[[1]] <- rbind(pointsMatrix[[1]], map[polyIndsRegion.i[k]][[1]])
          }
        else 
          {
          pointsMatrix <- map[polyIndsRegion.i]
          }
        pointsRegion.i <- sapply(X=pointsMatrix,
                                 FUN=pointsMatrixToPointsCharVector)

        ## j is the column region number
        for (j in (i + 1):nRegions)
        {
            ## which polygons belong to region number j?
            polyIndsRegion.j <- which(names(map) == regions[j])
            
            ## make a vector of all points (x, y) in the format "x&y" which belong to this region
            if(length(polyIndsRegion.j)>1)
              {
              pointsMatrix <- map[polyIndsRegion.j[1]]
              for(k in 2:length(polyIndsRegion.j))
                pointsMatrix[[1]] <- rbind(pointsMatrix[[1]], map[polyIndsRegion.j[k]][[1]])
              }
            else 
              {
              pointsMatrix <- map[polyIndsRegion.j]
              }
            pointsRegion.j <- sapply(X=pointsMatrix,
                                     FUN=pointsMatrixToPointsCharVector)

            ## now decide if region i and j share at least 2 common points in their polygons
            if (sum(pointsRegion.i %in% pointsRegion.j) >= npoints)
            {
                ## then they are neighbors!
                retMatrix[i, j] <-
                    retMatrix[j, i] <- - 1L
            }

            ## increment counter
            nPairsProcessed <- nPairsProcessed + 1L
            
            ## echo progress?
            if (nPairsProcessed %in% progressIterations)
            {
                ## echo percentage of processed pairs
                percentageProcessed <- floor((nPairsProcessed * 100L) / nPairs)
                cat(paste("progress: ", percentageProcessed, "%\n",
                          sep = ""))
            }            
        }
    }

    cat("Neighbor search finished.\n")

    ## add is.in relations
    surrounding <- attr(map, "surrounding")
    whichPolygonsAreInner <- which(sapply(surrounding, length) > 0L)   
    for(innerInd in whichPolygonsAreInner)
    {
        innerRegion <- names(map)[innerInd]
        outerRegion <- surrounding[[innerInd]]
        
        retMatrix[innerRegion, outerRegion] <-
            retMatrix[outerRegion, innerRegion] <- - 1L
    }

    ## on the diagonal, there are the negative row sums
    diag(retMatrix) <- - rowSums(retMatrix)

    ## finally return the matrix as a graph object
    class(retMatrix) <- "gra"
    return(retMatrix)
}


write.bnd <- function(map, file, replace=FALSE)
{
    if(! inherits(map,"bnd"))
        stop("argument 'map' is not an object of class 'bnd'")

    ## coercions, checks
    replace <- as.logical(replace)
    file <- as.character(file)
    stopifnot(identical(length(file), 1L),
              identical(length(replace), 1L))    

    ## check whether the file exists
    if(file.exists(file))
    {
        if(replace)
        {
            removeSucceeded <- file.remove(file)
            if(! removeSucceeded)
            {
                stop("file exists, but could not be removed")
            }
        } else {
            stop("specified file already exists")
        }        
    }

    myQuote <- function(string)
    {
        return(paste("\"", string, "\"",
                     sep=""))
    }
    
    ## names of the belonging regions
    belongingRegions <- names(map)
    
    ## no. of polygons
    nPolygons <- length(map)

    ## the surrounding list
    surrounding <- attr(map, "surrounding")
    
    for(i in seq_len(nPolygons))
    {
        dat <- map[[i]]
        dat <- dat[complete.cases(dat), ]
        
        temp <- paste(myQuote(belongingRegions[i]),
                      nrow(dat),
                      sep=",")
        write(temp, file, append=TRUE)
        
        if(length(outerRegionName <- surrounding[[i]]))
        {
            con <- paste("is.in",
                         myQuote(outerRegionName),
                         sep=",")
            write(con, file, append=TRUE)
        }
        write.table(dat, file, append=TRUE,
                    col.names=FALSE, row.names=FALSE,
                    sep=",", quote=FALSE)
    }

    return(invisible())
}


write.gra <- function(map, file, replace=FALSE)
{
  if(!inherits(map, "gra"))
    stop("Argument 'map' is not an object of class 'gra'!")

  ## check whether the file exists
  if(replace & file.exists(file))
    test <- file.remove(file)
  if(!replace & file.exists(file))
    stop("Specified file already exists!")

  ## names of districts
  districts <- as.integer(rownames(map))

  ## no. of regions
  S <- length(districts)
  write(S, file)

  ## loop over the regions
  for(i in 1:S){
    ## write name of the district
    write(districts[i], file, append=TRUE)
        
    ## write no. of neighbors
    write(map[i,i], file, append=TRUE)

    ## derive and write neighbors
    ind <- which(map[i, ] == -1) - 1
    write(ind, file, ncolumns = length(ind), append = TRUE)
  }

  return(invisible(NULL))
}


read.bnd <- function(file, sorted=FALSE)
{
    ## Liste wird erstellt: 1. Element enthält die 1. Spalte aus BND-Datei, 2. Element die 2. Spalte 
    data.raw <- scan(file,
                     what = list("", ""),
                     sep = ",",
                     quote = "")

    ## disable warnings and surely revert to old settings afterwards
    oldOptions <- options(warn = -1)
    on.exit(options(oldOptions))    

    ## Ursache für Warnungen: NAs werden bei jedem Regionsnamen und bei is.in erzeugt
    data.numeric <- lapply(data.raw,
                           as.numeric)

    ## revert now to old settings to show other warnings
    options(oldOptions)

    ## helper function
    unquote <- function(string)
    {
        return(gsub(pattern="\"",
                    replacement="",
                    x=string))
    }
    
    ## where do we have is.in's?
    whereIsIn <- which(data.raw[[1]] == "is.in") 

    ## so the surrounding names are
    surroundingNames <- unquote(data.raw[[2]][whereIsIn])

    ## so where do we have the region names of the polygons?
    whereRegionNames <- setdiff(which(is.na(data.numeric[[1]])),
                                whereIsIn)

    ## and we know how many polygons there are now
    nPolygons <- length(whereRegionNames)
    cat("Note: map consists of", nPolygons, "polygons\n")

    ## extract the region names
    belongingRegions <- unquote(data.raw[[1]][whereRegionNames])

    ## so we can already setup the regions vector
    regions <- unique(belongingRegions)
    cat("Note: map consists of", length(regions), "regions\n")
    
    ## extract the lengths of the polygons
    polyLengths <- data.numeric[[2]][whereRegionNames]

    ## find out which polygon number the is.in information belongs to
    enclosedPolygonsInds <- findInterval(x=whereIsIn,
                                         vec=whereRegionNames)

    ## so we can already setup the surrounding list
    surrounding <- replicate(n=nPolygons,
                             expr=character())
    for(i in seq_along(enclosedPolygonsInds))
    {
        surrounding[[enclosedPolygonsInds[i]]] <- surroundingNames[i]
    }
    
    ## now we can cbind and delete all not-point-stuff from the numeric data
    data.numeric <- cbind(data.numeric[[1]],
                          data.numeric[[2]])
    data.numeric <- na.omit(data.numeric)
    
    ## start processing the single polygons
    cat("Reading map ...")

    ## set up the list with the correct names
    map <- vector(mode="list", length=nPolygons)
    names(map) <- belongingRegions

    ## these are the start indices for the polygons
    ## (plus one index past the last polygon)
    startInds <- cumsum(c(1, polyLengths))
    
    ## so filling the list with the points is easy now
    for(k in seq_along(map)) {
        map[[k]] <- data.numeric[startInds[k]:(startInds[k+1] - 1), ]
        if(sum(map[[k]][1,] == map[[k]][polyLengths[k],]) != 2)
           warning(paste("Note: First and last point of polygon ",k," (region ",names(map)[k],") are not identical", sep=""), call. = FALSE)
    }

    ## processing finished
    cat(" finished\n")
    rm(data.numeric)
    
    ## sort the map?
    if(sorted){
        ## try conversion of region names to numbers
        numericNames <- as.numeric(names(map))

        ## decide which ordering to apply:
        ## are there some characters?
        newOrder <- 
            if(any(is.na(numericNames)))
            {
                cat("Note: regions sorted by name\n")
                order(names(map))
                
            } else {
                cat("Note: regions sorted by number\n")
                order(numericNames)                
            }

        ## order everything 
        map <- map[newOrder]
        surrounding <- surrounding[newOrder]
    }

    ## Bestimmung des Höhe-Breiten-Verhältnisses
    minima <- sapply(map, function(x){apply(x,2,min)})
    maxima <- sapply(map, function(x){apply(x,2,max)})
    
    minimum <- apply(minima,1,min)
    maximum <- apply(maxima,1,max)
    
    x.range <- maximum[1] - minimum[1]
    y.range <- maximum[2] - minimum[2]

    ## return the bnd object
    rval <- structure(map, class = "bnd",
      surrounding = surrounding, regions = regions)
    attr(rval, "asp") <- (y.range / x.range) / cos((mean(c(maximum[2] - minimum[2])) * pi) / 180)

    rval
}


read.gra <- function(file, sorted = FALSE, sep = " ")
{    
  ## read information from the graph file
  gfile <- strsplit(readLines(file), sep)

  ## get region names and neighbors
  if(length(gfile[[2]]) < 2) {
    ids <- rep(1:3, length = length(gfile) - 1)
    regions <- unlist(gfile[c(2:length(gfile))[ids == 1]])
    neighbors <- gfile[c(2:length(gfile))[ids == 3]]
    foo <- function(x) as.integer(x) + 1
  } else {
    regions <- sapply(gfile[-1], function(x) x[1])
    neighbors <- sapply(gfile[-1], function(x) x[-c(1:2)])
    foo <- function(x) as.integer(x)
  }
  neighbors <- lapply(neighbors, foo)

  ## extract the number of regions
  n <- length(regions)
  cat("Note: map consists of", n, "regions\n")

  cat("Creating adjacency matrix ...\n")

  ## create the adjacency matrix
  adjmat <- matrix(0, n, n)
  
  ## number of neighbors for each region
  diag(adjmat) <- sapply(neighbors, length)

  ## loop over the regions
  for(i in 1:n) {
    ## off-diagonal entries in adjmat
    ## Note: the indices are zero-based!
    if(length(neighbors[[i]]))
      adjmat[i, neighbors[[i]]] <- -1
  }

  colnames(adjmat) <- rownames(adjmat) <- regions

  cat("finished\n")
    
  if(sorted) {
    ## sort regions and adjacency matrix
    if(sum(is.na(as.numeric(regions))) == 0){
      regions <- as.numeric(regions)
      cat("Note: regions sorted by number\n") 
    } else cat("Note: regions sorted by name\n")
    ord <- order(regions)
    regions <- sort(regions)
    adjmat <- adjmat[ord, ord]
  }

  rownames(adjmat) <- colnames(adjmat) <- regions

  class(adjmat) <- c("gra", "matrix")
  return(adjmat)
}


shp2bnd <- function(shpname, regionnames, check.is.in = TRUE)
{   
  require("shapefiles")

    ## safe coercions ...
    shpname <- as.character(shpname)
    regionnames <- as.character(regionnames)
    check.is.in <- as.logical(check.is.in)

    ## ... and checks
    stopifnot(identical(length(shpname), 1L),
              length(regionnames) >= 1L,
              identical(length(check.is.in), 1L))
    
    ## now read the shapefile information
    shp <- shapefiles::read.shapefile(shpname)
    dbf <- shapefiles::read.dbf(paste(shpname,".dbf",sep=""))

    ## extract names of the regions:
    regionnames <-
        ## if there is only one string, then we assume it is the variable name
        ## in the database of the shape file
        if(identical(length(regionnames), 1L))
        {
            ## so get that variable
            as.character(dbf$dbf[regionnames][[1L]])
        } else {
            ## we stay at the given names
            regionnames
        }
    
    ## delete commas in names of the regions
    newRegionnames <- gsub(pattern=",",
                           replacement="",
                           x=regionnames,
                           fixed=TRUE)

    ## check if commas have been deleted and if so, issue a warning
    if(! identical(regionnames, newRegionnames))
        warning(simpleWarning("commas in names of the regions have been deleted"))

    ## overwrite old region names with new region names
    regionnames <- newRegionnames
    rm(newRegionnames)

    ## split data into closed polygons.
    ## we need:
    
    ## storage for the new names and the new polygons
    polyList <- list()
    ## and the corresponding original indexes
    originalRegionNumber <- integer()

    cat("Reading map ...")
    
    ## now process all original regions
    for(i in seq_along(regionnames))
    {
        ## this is temporary storage (originally a data.frame with X and Y columns):
        temppoly <- as.matrix(shp$shp$shp[[i]]$points)
        dimnames(temppoly) <- NULL

        ## as long as there are still points to be processed
        while((nPoints <- nrow(temppoly)) > 0L)
        {
            ## where does the first point occur in the data at the second time?
            endIndex <- which((temppoly[- 1L, 1] == temppoly[1L, 1]) &
                              (temppoly[- 1L, 2] == temppoly[1L, 2])) + 1L

            ## take the first next occurrence, or the last point if the polygon is not closed
            endIndex <- 
                if(length(endIndex) > 0L)
                {
                    endIndex[1L]
                } else {
                    nPoints
                }           

            ## the range of this polygon
            polyRange <- 1L:endIndex

            ## this was index i
            originalRegionNumber <- c(originalRegionNumber,
                                      i)
            
            ## save the polygon
            polyList <- c(polyList,
                          list(temppoly[polyRange, ]))
            ## list is necessary so that c(list(), data.frame(..)) is a one-element list,
            ## and not a list with the variables of the data.frame as elements

            ## and delete this part from temporary storage
            temppoly <- temppoly[- polyRange, ]
        }
    }

    cat(" finished\n")
    
    ## so how many polygons do we have now?
    nPolys <- length(polyList)

    cat("Note: map consists originally of", nPolys, "polygons\n")

    ## here is the parallel list of the surrounding region names of single polygons
    surrounding <- replicate(n=nPolys,
                             expr=character()) ## until now no region names anywhere!
    
    ## check for polygons contained in another polygon?
    if(check.is.in)
    {      
        ## get dimensions of all polygons
        dims <- sapply(polyList, nrow)

        ## save here which polygons should be removed, because they are boundaries
        ## to polygons lying inside
        rmcheck <- logical(nPolys)

        ## save here the indexes of the polygons which have already been matched/processed.
        ## these must not be processed again!
        whichWereProcessed <- integer()

        ## process each polygon i
        for(i in seq_len(nPolys))
        {
            ## if we had processed this already
            if(i %in% whichWereProcessed)
            {
                ## go on to the next polygon
                next
            } else {
                ## add i to processed ones
                whichWereProcessed <- union(whichWereProcessed,
                                            i)
            }
            
            ## which polygons have same number of points as the current?
            sameDimsInds <- setdiff(which(dims == dims[i]),
                                    whichWereProcessed) ## but without the already processed ones     

            ## process all polygons j with same dims as polygon i
            ## (this works as a hash)
            for(j in sameDimsInds)
            {
                ## compute squared distance of polygon_i and reversed polygon_j
                reverseInds <- dims[i]:1L
                squaredDistance <- sum( (polyList[[i]] -
                                         polyList[[j]][reverseInds, ])^2 )

                ## if it is small enough
                if(squaredDistance < 1e-5)
                {                       
                    ## find out which is the outer one
                    outer <- inner <- 0L
                    if(.ringDirxy(polyList[[j]]) < 0)
                    {
                        outer <- j
                        inner <- i
                    } else {
                        outer <- i
                        inner <- j
                    }

                    ## remove the outer polygon
                    rmcheck[outer] <- TRUE

                    ## and add the information in which region it is lying
                    ## (each polygon can only lie in 1 other region, of course...)
                    surrounding[[inner]] <- regionnames[originalRegionNumber[outer]]
                }

                ## we have processed j
                whichWereProcessed <- union(whichWereProcessed,
                                            j)
            }            
        }

        ## we have processed all polygons, and can remove the unnecessary ones
        polyList <- polyList[! rmcheck]
        originalRegionNumber <- originalRegionNumber[! rmcheck]
        surrounding <- surrounding[! rmcheck]

        cat("Note: After removing unnecessary surrounding polygons, the map consists of",
            length(polyList), "polygons\n")          
    }    

    ## add the original region names to the polygons list as names
    names(polyList) <- regionnames[originalRegionNumber]

    ## the new unique regions
    regions <- unique(names(polyList))
    cat("Note: map consists of", length(regions), "regions\n")
    
    ## compute relation of height to width (for plotting etc)
    minima <- sapply(polyList, function(x){apply(x,2,min)})
    maxima <- sapply(polyList, function(x){apply(x,2,max)})
    
    minimum <- apply(minima,1,min)
    maximum <- apply(maxima,1,max)
    
    x.range <- maximum[1] - minimum[1]
    y.range <- maximum[2] - minimum[2]
    
    height2width <- round(y.range / x.range, digits=2)
    
    ## now return the bnd object
    return(structure(polyList,
                     class="bnd",
                     height2width=height2width,
                     surrounding=surrounding,
                     regions=regions))
}



nb2gra <- function(nbObject)
{
    ## check if S3 class of nbObject is "nb"
    stopifnot(inherits(x=nbObject,
                       what="nb"))

    ## convert to (negative) binary neighbors matrix
    regionIds <- attr(nbObject, "region.id")
    ret <- matrix(data=0,
                  nrow=length(regionIds),
                  ncol=length(regionIds),
                  dimnames=
                  list(regionIds,
                       regionIds))

    for(i in seq_along(nbObject))
    {
        ret[i, nbObject[[i]]] <- - 1
    }


    ## and to gra format
    diag(ret) <- - rowSums(ret)
    class(ret) <- "gra"

    return(ret)
}


gra2nb <- function(graObject)
{
    ## check if S3 class of nbObject is "gra"
    stopifnot(inherits(x=graObject,
                       what="gra"))

    ## save region names and delete them
    ## (so that the list below will not have names attached)
    regionNames <- rownames(graObject)
    dimnames(graObject) <- NULL
    
    ## make list of neighbors
    ret <- apply(graObject,
                 MARGIN=1,
                 FUN=function(row) which(row == -1))

    ## attach necessary attributes
    ret <- structure(ret,
                     class="nb",
                     region.id=regionNames,
                     call=match.call(),
                     type="queen",
                     sym=TRUE)

    ## and return the nb object
    return(ret)
}



bnd2sp <- function(bndObject)
{
  require("sp")

    ## check if S3 class of bndObject is "bnd"
    stopifnot(inherits(x=bndObject,
                       what="bnd"))

    ## extracts
    bndNames <- names(bndObject)
    regions <- unique(bndNames)
    bndAttributes <- attributes(bndObject)
    
    ## close all polygons (last coordinates must match first coordinates)
    bndObject <- lapply(bndObject,
                        FUN=function(polygon){
                            if(! isTRUE(identical(polygon[1, ],
                                                  polygon[nrow(polygon), ])))
                            {
                                rbind(polygon,
                                      polygon[1, ])
                            } else {
                                polygon
                            }
                        })
    
    ## set up return list
    ret <- list()
    
    ## process all unique regions
    for(id in regions)
    {
        ## which polygons belong to this region?
        idMatches <- which(id == bndNames)

        ## convert these polygons to Polygon class objects
        idPolygons <- lapply(bndObject[idMatches],
                             FUN=sp::Polygon,
                             hole=FALSE)

        ## add the Polygons object with these Polygon parts to return list
        ret[[id]] <- sp::Polygons(srl=idPolygons,
                                  ID=id)        
    }

    ## add holes of inner polygons to outer regions
    surrounding <- bndAttributes$surrounding
    whichAreInner <- which(sapply(surrounding, length) > 0L)
    for(innerInd in whichAreInner)
    {
        ## get the hole
        hole <- sp::Polygon(coords=bndObject[[innerInd]],
                            hole=TRUE)

        ## get outer polys list
        outerId <- surrounding[[innerInd]]
        outerPolys <- ret[[outerId]]@Polygons

        ## add the hole to outer polys
        outerPolys <- c(outerPolys,
                        hole)

        ## write back extended outer polys list as new Polygons object with same ID as before 
        ret[[outerId]] <- sp::Polygons(srl=outerPolys,
                                       ID=outerId)
    }
    
    ## convert list of Polygons to a SpatialPolygons object and return that
    ret <- sp::SpatialPolygons(Srl=ret)
    return(ret)
}


sp2bnd <- function(spObject,
  ## object of class SpatialPolygons (or specializations)
  regionNames=sapply(spObject@polygons, slot, "ID"),
  ## character vector of region names
  ## (parallel to the Polygons list in spObject)
  height2width=round(diff(sp::bbox(spObject)[2, ]) / diff(sp::bbox(spObject)[1, ]), 2),
  ## ratio of height to width
  epsilon=sqrt(.Machine$double.eps))
  ## how much can two polygons differ (in maximumn
  ## squared Euclidean distance) and still match each other?
    
{
  require("sp")

  ## check if S4 class of spObject is "SpatialPolygons"
  stopifnot(is(object = spObject, class2 = "SpatialPolygons"))

  ## extracts
  spObject <- sp::polygons(spObject)  ## now surely a SpatialPolygons object
  spList <- spObject@polygons         ## discard other slots
  nRegions <- length(spList)

  ## check if number of regions matches with the length of regionNames etc
  stopifnot(is.character(regionNames),
    identical(length(regionNames), nRegions),
    height2width > 0)

  ## set up return and holes list
  ret <- list()
  holes <- list()    

  ## number of polygons and holes already processed
  numPolysProcessed <- 0
  numHolesProcessed <- 0
    
  ## process each region
  for(regionIterator in seq_along(spList)) {
    thisRegion <- spList[[regionIterator]]@Polygons
        
        ## process each Polygon in this region
        for(polygonObject in thisRegion)
        {
            ## if it is a hole, put it in holes, else in ret.
            ## the name is set to the region name so we know from which region this
            ## polygon stems.
            if(polygonObject@hole)
            {
                ## increment hole counter
                numHolesProcessed <- numHolesProcessed + 1

                ## and correct invariant
                holes[[numHolesProcessed]] <- sp::coordinates(polygonObject)
                names(holes)[numHolesProcessed] <- regionNames[regionIterator]
            } else {
                ## increment Polygon counter
                numPolysProcessed <- numPolysProcessed + 1

                ## and correct invariant
                ret[[numPolysProcessed]] <- sp::coordinates(polygonObject)
                names(ret)[numPolysProcessed] <- regionNames[regionIterator]
            }
        }
    }
    ## sanity check
    stopifnot(all.equal(length(ret), numPolysProcessed),
              all.equal(length(holes), numHolesProcessed))

    ## now process all holes:

    ## set up surrounding list
    surrounding <- replicate(n=numPolysProcessed,
                             character())
    
    ## use number of coordinates as hash for quicker search for the embedded region
    polyDims <- sapply(ret, nrow)
    holeDims <- sapply(holes, nrow)

    for(i in seq_along(holes))
    {
        ## hash lookup
        potentialMatchesInds <- which(holeDims[i] == polyDims)

        ## now more precise search in these potential matches
        matchFound <- FALSE
        for(j in potentialMatchesInds)
        {
            ## decide
            thisHole <- holes[[i]]
            thisRegion <- ret[[j]]

            squaredEuclideanDistances <-
                rowSums((thisHole[rev(seq_len(nrow(thisHole))), ] - thisRegion)^2) 
            doesMatch <- max(squaredEuclideanDistances) < epsilon

            ## if it matches, update the surrounding data
            ## and break out of the for loop
            if(doesMatch)
            {
                matchFound <- TRUE

                surrounding[[j]] <- names(holes)[i]
                
                ## we can proceed with the next hole:
                break
            }
        }

        ## echo a warning if a hole has no match
        if(! matchFound)
        {
            warning(simpleWarning(paste("No match found for hole in region",
                                        names(holes)[i])))
        }
    }

    ## finally collect information and return the bnd object
    ret <- structure(ret,
                     surrounding=surrounding,
                     height2width=height2width,
                     class="bnd")
    return(ret)
}


add.neighbor <- function(map, region1, region2)
{
   if(! inherits(map,"gra"))
      stop("Argument 'map' is not an object of class 'gra'!")

   # names of districts
   districts <- rownames(map)

   # find index of the regions with changes
   ind1 <- which(districts == region1)
   ind2 <- which(districts == region2)
  
   # modify neighborhoopd structure
   map[ind1,ind2] <- map[ind2,ind1] <- -1
  
   # and adjust no. of neighbors
   map[ind1,ind1] <- -sum(map[ind1, -ind1])
   map[ind2,ind2] <- -sum(map[ind2, -ind2])

   class(map) <- "gra"
   return(map)
}
