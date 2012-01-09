run.bayesx <-
function(prg = NULL, verbose = TRUE, ...)
{
  os.type <- .Platform$OS.type == "windows"
  if(os.type) {
    bin <- shQuote(system.file(package = "BayesXsrc", "libs", .Platform$r_arch, "BayesX.exe"))
  } else {
    bin <- shQuote(system.file(package = "BayesXsrc", "libs", .Platform$r_arch, "BayesX"))
  }
  if(is.null(prg)) {
    if(!os.type)
      bin <- paste(shQuote(file.path(R.home(),"bin","R")), "CMD", bin)
    system(bin, ...)
    return(invisible(NULL))
  } else {
    fileo <- getwd()
    prg <- path.expand(prg)
    dir <- dirname(prg)
    prg.name <- basename(prg)
    setwd(dir)
    output <- file.exists("output")
    temp <- file.exists("temp")
    win <- FALSE
    ptm <- proc.time()
    if(!verbose) {
      if(.Platform$OS.type == "unix") {
        log <- try(system(paste(bin, " ", dir, "/", prg.name, 
          " > ", dir, "/bayesx.log", sep = ""), intern = FALSE, ...))
      } else {
        log <- try(system(paste(bin, " ", dir, "/", prg.name, sep = ""), 
          intern = TRUE, show.output.on.console = FALSE, ...))
        win <- TRUE
      }
    } else log <- try(system(paste(bin, " ", dir, "/", prg.name, sep = ""), intern = FALSE, ...))
	  if(log != 0 && !win)
      warning("problem processing BayesX!")
    if(!verbose && .Platform$OS.type == "unix")	
      log <- readLines(paste(dir, "/bayesx.log", sep = ""))
    now <- proc.time()
    runtime <- now - ptm
    runtime <- runtime
    if(verbose)
      cat("Total run time was:", runtime[3L], "sec\n")
    if(!output)
      output <- file.exists("output")
    if(output)
      try(suppressWarnings(file.remove("output")), silent = TRUE)
    if(!output)
      output <- file.exists("temp")
    if(output)
      try(suppressWarnings(file.remove("temp")), silent = TRUE)
    try(setwd(fileo), silent = TRUE)
    logerror <- FALSE
    for(logline in log)
      if(grepl("ERROR", logline))
        logerror <- TRUE
    if(logerror)
      warning("an error occurred during runtime of BayesX, please check the BayesX logfile!")
    return(invisible(list(log = log, runtime = runtime)))
  }
}
