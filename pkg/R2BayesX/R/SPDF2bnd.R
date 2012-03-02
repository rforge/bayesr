SPDF2bnd <- function(map)
{
  if(!is(map, "SpatialPolygonsDataFrame"))
    stop("map must be of class 'SpatialPolygonsDataFrame'!")
  rval <- list(); poly.names <- NULL
  for(i in 1:length(slot(map, "polygons"))) {
    rval[[i]] <- slot(slot(slot(map, "polygons")[[i]], "Polygons")[[1]], "coords")
    poly.names <- c(poly.names, slot(slot(map, "polygons")[[i]], "ID"))
  }
  names(rval) <- poly.names
  class(rval) <- c("bnd", "list")
  bb <- bbox(map)
  ## attr(rval, "asp") <- as.numeric(mapasp(map, bb[1L,], bb[2L,]))
  rval
}
