s4class <-
function(x)
{
  if(grepl("_random.", x, fixed = TRUE))
    cx <- "random.bayesx"
  if(grepl("_pspline.", x, fixed = TRUE))
    cx <- "sm.bayesx"
  if(grepl("_season.", x, fixed = TRUE))
    cx <- "sm.bayesx"
  if(grepl("_rw.", x, fixed = TRUE))
    cx <- "sm.bayesx"
  if(grepl("_spatial.", x, fixed = TRUE))
    cx <- "mrf.bayesx"
  if(grepl("_geospline.", x, fixed = TRUE))
    cx <- "geo.bayesx"
  if(grepl("_geokriging.", x, fixed = TRUE))
    cx <- "geo.bayesx"
  if(grepl("_logbaseline.", x, fixed = TRUE))
    cx <- "sm.bayesx"
  if(grepl("_kriging.", x, fixed = TRUE))
    cx <- "sm.bayesx"

  return(cx)
}

s4bs <-
function(x)
{
  if(grepl("_random", x))
    bs <- "re"
  if(grepl("_pspline", x))
    bs <- "ps"
  if(grepl("_season", x))
    bs <- "season"
  if(grepl("_rw", x))
    bs <- "rw"
  if(grepl("_spatial", x))
    bs <- "mrf"
  if(grepl("_geospline", x))
    bs <- "gs"
  if(grepl("_geokriging", x))
    bs <- "gk"
  if(grepl("_logbaseline", x))
    bs <- "bl"
  if(grepl("_kriging", x))
    bs <- "kr"

  return(bs)
}
