s4class <-
function(x)
{
  if(grepl("_random", x))
    cx <- "random.bayesx"
  if(grepl("_pspline", x))
    cx <- "sm.bayesx"
  if(grepl("_season", x))
    cx <- "sm.bayesx"
  if(grepl("_season", x))
    cx <- "sm.bayesx"
  if(grepl("_rw", x))
    cx <- "sm.bayesx"
  if(grepl("_spatial", x))
    cx <- "mrf.bayesx"
  if(grepl("_geospline", x))
    cx <- "geo.bayesx"
  if(grepl("_geokriging", x))
    cx <- "geo.bayesx"
  if(grepl("_logbaseline", x))
    cx <- "sm.bayesx"
  if(grepl("_kriging", x))
    cx <- "sm.bayesx"

  return(cx)
}

