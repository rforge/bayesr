library("zoo")

## Paper: http://arxiv.org/pdf/1503.07709.pdf
## Data: http://data.london.gov.uk/dataset/london-fire-brigade-incident-records
## Note: renamed data set, saved as .csv file...
d <- read.csv("~/data/LondonFire/data2013-2016.csv",
  header = TRUE, stringsAsFactors = FALSE)
d <- subset(d, PropertyCategory == "Dwelling" & IncidentGroup == "Fire")
d$DateOfCall <- as.Date(d$DateOfCall, "%d-%b-%y")

LondonFire <- list()
LondonFire$year <- as.integer(format(as.yearmon(d$DateOfCall), "%Y"))
LondonFire$month <- as.integer(format(as.yearmon(d$DateOfCall), "%m"))
LondonFire$day <- as.integer(format(as.yearmon(d$DateOfCall), "%d"))
time <- as.POSIXlt(d$TimeOfCall, format = "%H:%M:%S")
LondonFire$daytime <- time$hour + time$min / 60 + time$sec / (60^2)
LondonFire$east <- as.numeric(d$Easting_rounded)
LondonFire$north <- as.numeric(d$Northing_rounded)
LondonFire$arrivaltime <- as.numeric(d$FirstPumpArriving_AttendanceTime) / 60
LondonFire$id <- as.factor(d$IncidentNumber)
LondonFire <- na.omit(as.data.frame(LondonFire))

LondonFire <- subset(LondonFire, year == 2015)


## Extract fire stations address from:
## http://www.london-fire.gov.uk/A-ZFireStations.asp
## Copy & paste to spreadsheet and save as .csv file.
fs <- readLines("~/data/LondonFire/FireStations.csv")
fs <- gsub('\"', "", fs, fixed = TRUE)
fs <- fs[!grepl("redeveloped", fs)]
i <- grep("Borough", fs)
sn <- as.character(na.omit(c(fs[1], fs[i + 1])))
sa <- NULL
for(i in 1:(length(sn) - 1)) {
  j1 <- which(fs == sn[i]) + 1
  j2 <- which(fs == sn[i + 1]) - 2
  sa <- c(sa, paste(c("London", rev(fs[j1:j2])), collapse = ","))
}
j <- which(fs == sn[length(sn)]) + 1
sa <- c(sa, paste(c("London", rev(fs[j:(length(fs) - 1)])), collapse = ","))

trim <- function(x) {
  gsub("^\\s+|\\s+$", "", x)
}

sa <- trim(sa)
sn <- trim(sn)

LondonFStations <- data.frame("name" = sn, "address" = sa, stringsAsFactors = FALSE)


## Add fire stations coordinates.
library("RCurl")
library("rjson")

url <- function(address, return.call = "json", sensor = "false") {
  root <- "http://maps.google.com/maps/api/geocode/"
  u <- paste(root, return.call, "?address=", address, "&sensor=", sensor, sep = "")
  return(URLencode(u))
}


geoCode <- function(address,verbose=FALSE) {
  co <- NULL
  for(j in seq_along(address)) {
    cat("searching address", j, "of", length(address), "\n")
    if(verbose) cat(address[j],"\n")
    u <- url(address[j])
    doc <- getURL(u)
    x <- fromJSON(doc)
    if(x$status=="OK") {
      lat <- x$results[[1]]$geometry$location$lat
      lng <- x$results[[1]]$geometry$location$lng
      location_type  <- x$results[[1]]$geometry$location_type
      formatted_address  <- x$results[[1]]$formatted_address
      co <- rbind(co, data.frame("lon" = lng, "lat" = lat,
        "loc_type" = location_type, "address" = formatted_address,
        stringsAsFactors = FALSE))
      Sys.sleep(runif(1, 0.2, 0.6))
    } else {
      warning(paste("could not find location ", address[j], "!", sep = ""))
      cbind(co, data.frame("lon" = NA, "lat" = NA, "loc_type" = NA, "address" = NA))
    }
  }
  return(co)
}


## Get the stations lon,lat.
library("sp")
co <- geoCode(LondonFStations$address)
LondonFStations <- cbind(LondonFStations, co[, c("lon", "lat")])
coordinates(LondonFStations) <- c("lon", "lat")
proj4string(LondonFStations) <- CRS("+init=epsg:4326")


## Transform coordinates.
xy <- LondonFire[, c("east", "north")]
coordinates(xy) <- c("east", "north")
proj4string(xy) <- CRS("+init=epsg:27700")
res <- spTransform(xy, CRS("+init=epsg:4326"))
res <- as.data.frame(res)
names(res) <- c("lon", "lat")
LondonFire$lon <- res$lon
LondonFire$lat <- res$lat


## Calculate fire stations intensity.
library("spatstat")
library("raster")
chull <- convexhull.xy(as.matrix(rbind(LondonFire[ , c("lon", "lat")],
  coordinates(LondonFStations))))
win <- expandwinPerfect(chull, TRUE, 1)
fsppp <- ppp(x = LondonFStations$lon, LondonFStations$lat, window = win)
fsintens <- density.ppp(fsppp, bw.diggle)
fsintens <- raster(fsintens)
proj4string(fsintens) <- CRS("+init=epsg:4326")
LondonFire$fsintens <- extract(fsintens, as.matrix(LondonFire[ , c("lon", "lat")]))


## Finalize.
LondonFire <- LondonFire[, c("arrivaltime", "daytime", "fsintens", "lon", "lat")]
coordinates(LondonFire) <- c("lon", "lat")
proj4string(LondonFire) <- CRS("+init=epsg:4326")


## Get London boroughs.
## http://data.london.gov.uk/dataset/statistical-gis-boundary-files-london
London.maps <- function()
{
  require("maptools")

  owd <- getwd()
  dir.create(tdir <- tempfile())

  download.file("https://files.datapress.com/london/dataset/statistical-gis-boundary-files-london/statistical-gis-boundaries-london.zip", file.path(tdir, "London.zip"))

  setwd(tdir)
  unzip(file.path(tdir, "London.zip"), exdir = tdir)

  lb <- readShapePoly(file.path(tdir, "statistical-gis-boundaries-london",
    "ESRI", "London_Borough_Excluding_MHW.shp"), proj4string = CRS("+init=epsg:27700"))
  lb <- spTransform(lb, CRS("+init=epsg:4326"))

  lb0 <- readShapePoly(file.path(tdir, "statistical-gis-boundaries-london",
    "ESRI", "London_Ward_CityMerged.shp"), proj4string = CRS("+init=epsg:27700"))
  lb0 <- spTransform(lb0, CRS("+init=epsg:4326"))
  lb0 <- unionSpatialPolygons(lb0, rep(1L, length = length(lb0)))

  lb <- as(lb, "SpatialPolygons")
  lb0 <- as(lb0, "SpatialPolygons")

  setwd(owd)

  return(list("Boroughs" = lb, "Boundaries" = lb0))
}

London <- London.maps()
LondonBoroughs <- London$Boroughs
LondonBoundaries <- London$Boundaries


## Save the data.
library("bamlss")

save_data(LondonFire, LondonFStations, LondonBoroughs, LondonBoundaries,
  file = "~/svn/bayesr/pkg/bamlss/data/LondonFire.rda")

