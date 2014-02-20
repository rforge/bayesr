## Function to save data in .rda format with
## minimum storing size.
save_data <- function(..., file, envir = parent.frame())
{
  dir.create(tdir <- tempfile())
  on.exit(unlink(tdir))
  bf <- basename(file)
 
  compress <- c("gzip", "bzip2", "xz")
  size <- NULL
  for(j in compress) {
    tf <- file.path(tdir, paste(j, bf, sep = "-"))
    save(..., file = tf, compress = j, envir = envir)
    size <- c(size, file.info(tf)$size)
  }

  print(data.frame("compress" = compress, "size" = size))

  compress <- compress[which.min(size)]
  save(..., file = file, compress = compress, envir = envir)
}


## Munich rent index.
data_MunichRent <- function(dir = NULL)
{
  if(is.null(dir))
    dir <- "~/svn/bayesr/pkg/BayesR/data"
  dir <- path.expand(dir)

  dpath <- "http://www.stat.uni-muenchen.de/~kneib/regressionsbuch/download/mietspiegel99.raw"
  dat <- read.table(dpath, header = TRUE)
  rent99 <- list()
  rent99$rent <- dat$miete
  rent99$rentsqm <- dat$mieteqm
  rent99$area <- dat$flaeche
  rent99$yearc <- dat$bjahr
  rent99$bath <- as.factor(dat$bad)
  levels(rent99$bath) <- c("standard", "premium")
  rent99$kitchen <- as.factor(dat$kueche)
  levels(rent99$kitchen) <- c("standard", "premium")
  rent99$district <- dat$bezv
  rent99$location <- as.factor(dat$lage)
  levels(rent99$location) <- c("average", "good", "top")
  rent99$cheating <- as.factor(dat$zh)
  levels(rent99$cheating) <- c("no", "yes")
  rent99 <- as.data.frame(rent99)
  rent99 <- rent99[order(rent99$district), ]

  nenv <- new.env()
  assign("rent99", rent99, envir = nenv)
  save_data(rent99, file = file.path(dir, "rent99.rda"), envir = nenv)

  dpath <- "http://www.stat.uni-muenchen.de/~kneib/regressionsbuch/download/muenchen.bnd"
  MunichBnd <- read.bnd(dpath)
  nm <- names(MunichBnd)
  MunichBnd <- MunichBnd[nm[order(as.integer(nm))]]
  attr(MunichBnd, "asp") <- 1.1

  assign("MunichBnd", MunichBnd, envir = nenv)
  save_data(MunichBnd, file = file.path(dir, "MunichBnd.rda"), envir = nenv)

  invisible(NULL)
}


## Patent opposition data.
data_Patent <- function(dir = NULL)
{
  if(is.null(dir))
    dir <- "~/svn/bayesr/pkg/BayesR/data"
  dir <- path.expand(dir)

  dpath <- "http://www.stat.uni-muenchen.de/~kneib/regressionsbuch/download/patentdata.raw"
  dat <- read.table(dpath, header = TRUE)
  patent <- list()
  patent$opposition <- factor(dat$opp, levels = c(0, 1), labels = c("no", "yes"))
  patent$biopharm <- factor(dat$biopharm, levels = c(0, 1), labels = c("no", "yes"))
  patent$USA2 <- factor(dat$ustwin, levels = c(0, 1), labels = c("no", "yes"))
  patent$holder <- factor(dat$patus, levels = c(0, 1), labels = c("EU", "USA"))
  patent$GSGB <- factor(dat$patgsgr, levels = c(0, 1), labels = c("no", "yes"))
  patent$year <- as.integer(dat$year)
  patent$ncitations <- as.integer(dat$ncit)
  patent$ncountry <- as.integer(dat$ncountry)
  patent$nclaims <- as.integer(dat$nclaims)

  nenv <- new.env()
  assign("patent", patent, envir = nenv)
  save_data(patent, file = file.path(dir, "patent.rda"), envir = nenv)

  invisible(NULL)
}

