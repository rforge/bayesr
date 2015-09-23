library("zoo")

## Paper: http://arxiv.org/pdf/1503.07709.pdf
## Data: http://data.london.gov.uk/dataset/london-fire-brigade-incident-records
## Note: renamed data sets...
d <- read.csv("~/data/LondonFire/data2012-2015.csv",
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

FireStations <- data.frame("name" = sn, "address" = sa, stringsAsFactors = FALSE)


## Add fire stations coordinates.
library("RCurl")
library("rjson")

key <- "AIzaSyBEtRUvHLe8Ah6N022HHXuQKRt5rDvsxGQ"

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


co <- geoCode(FireStations$address)

FireStations <- cbind(FireStations, co[, c("lon", "lat")])


## Transform coordinates.
library("sp")
xy <- LondonFire[, c("east", "north")]
coordinates(xy) <- c("east", "north")
proj4string(xy) <- CRS("+init=epsg:27700")
res <- spTransform(xy, CRS("+init=epsg:4326"))
res <- as.data.frame(res)
names(res) <- c("lon", "lat")
LondonFire$lon <- res$lon
LondonFire$lat <- res$lat


## Calculate fire stations intensity
library("spatstat")
library("raster")
chull <- convexhull.xy(as.matrix(rbind(LondonFire[ , c("lon", "lat")],
  FireStations[, c("lon", "lat")])))
win <- expandwinPerfect(chull, TRUE, 1)
fsppp <- ppp(x = FireStations$lon, FireStations$lat, window = win)
fsintens <- density.ppp(fsppp, bw.diggle)
fsintens <- raster(fsintens)
proj4string(fsintens) <- CRS("+init=epsg:4326")
LondonFire$fsintens <- extract(fsintens, as.matrix(LondonFire[ , c("lon", "lat")]))


## Finalize.
LondonFire <- LondonFire[, c("year", "month", "day", "daytime",
  "arrivaltime", "id", "lon", "lat", "fsintens")]

save(LondonFire, FireStations,
  file = "~/svn/bayesr/pkg/bamlss/data/LondonFire.rda",
  compress = "xz")

