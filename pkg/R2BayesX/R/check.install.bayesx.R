check.install.bayesx <-
function(bin = getOption("bayesx.bin"), verbose = FALSE)
{
  options(warn = -1)
  ## trivial test program
  temp.dir <- tempdir()
  bayesx.test.prg <- "exit\n"
  testfile <- paste(temp.dir, "/bayesx.test.prg", sep = "")
  writeLines(bayesx.test.prg,testfile)
  ok <- TRUE
  if(!is.null(bin))
    if(!file.exists(bin))
      ok <- FALSE

  ## if bin unspecified try R2BayesX/bin or direct call
  if(is.null(bin))  {
    base_name <- if(.Platform$OS.type == "unix") "BayesX" else "bayesx.exe"
    if(!file.exists(bin <- file.path(.find.package("R2BayesX"), "bin", base_name))) 
      bin <- base_name 
  }

  ## check
  if(ok) {
    owd <- getwd()
    setwd(temp.dir)
    try(suppressWarnings(dir.create(paste(temp.dir, "/output", sep = ""))), silent = TRUE)
    try(suppressWarnings(dir.create(paste(temp.dir, "/temp", sep = ""))), silent = TRUE)
    ok <- try(suppressWarnings(system(paste(bin, testfile), intern = TRUE)), silent = TRUE)
    if(length(ok) > 0 && !inherits(ok, "try-error")) {
      try(Sys.chmod("output", mode = "0777"), silent = TRUE)
      try(Sys.chmod("temp", mode = "0777"), silent = TRUE)
      try(suppressWarnings(file.remove("output")), silent = TRUE)
      try(suppressWarnings(file.remove("temp")), silent = TRUE)
    }
    if(!is.null(owd))
      setwd(owd)
  }
  if(length(ok) > 0) {
    if(ok == "> exit")
      ok <- TRUE
    else
      ok <- FALSE
  } else ok <- FALSE

  ## report results
  if(verbose) {
    if(ok) 
      cat(sprintf("BayesX found in: %s\n", bin))
    else
      warning("BayesX not found: please install BayesX or (if already installed) verify bin path!")
  }
  options(warn = 0)

  if(ok) 
    return(bin) 
  else 
    return(NULL)
}

