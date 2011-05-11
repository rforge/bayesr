run.bayesx <-
function(dir, prg.name = "bayesx.estim.input.prg", 
  verbose = FALSE, bin = getOption("bayesx.bin"))
{
  ok <- check.install.bayesx(bin, verbose = FALSE) 
  if(is.null(ok))
    stop("please install BayesX or verify bin path, also see function install.bayesx()")
  fileo <- getwd()
  setwd(dir)
  output <- file.exists("output")
  temp <- file.exists("temp")
  win <- FALSE
  ptm <- proc.time()
  if(!verbose) {
    if(.Platform$OS.type == "unix") {
      log <- try(system(paste(ok, " ", dir, "/", prg.name, 
        " > ", dir, "/bayesx.log", sep = ""), intern = FALSE))
    } else {
      log <- try(system(paste(ok, " ", dir, "/", prg.name, sep = ""), 
        intern = TRUE, show.output.on.console = FALSE))
      win <- TRUE
    }
  } else log <- try(system(paste(ok, " ", dir, "/", prg.name, sep = ""), intern = FALSE))
	if(log != 0 && !win)
    stop("BayesX could not be processed (maybe verify bin path)!")
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
  if(logerror) {
    warning(paste("an error occurred during runtime of BayesX, please check the logfile,", 
      "see the help site of function bayesx_logfile()!"))
  }

  return(invisible(list(log = log, runtime = runtime)))
}

