install.bayesx <- function(inst.dir = NULL, source.dir = NULL)
{
  options(warn = -1)
  unixos <- .Platform$OS.type == "unix"
  base_name <- if(unixos) "BayesX" else "bayesx.exe"
  if(!is.null(inst.dir))
    inst.dir <- path.expand(inst.dir)
  if(!is.null(source.dir))
    source.dir <- path.expand(source.dir)

  ## check if bin folder in package may be created
  lib.dir <- .find.package("R2BayesX")
  if(is.null(inst.dir))
    bin.dir <- paste(lib.dir, "/bin", sep = "")
  else
    bin.dir <- inst.dir
  ok <- file.exists(bin.dir)
  if(!ok && is.null(inst.dir))
    ok <- try(suppressWarnings(dir.create(bin.dir)), silent = TRUE)
  if(inherits(ok, "try-error"))
    ok <- FALSE

  ## if not need to supply inst.dir
  if(!ok) {
    if(is.null(inst.dir)) {
      cat("the install directory has not been supplied yet, please specify now!\n")
      bin.dir <- readline()
    }
    ok <- file.exists(bin.dir)
    if(!ok) {
      options(warn = 0)
      stop(paste("wrong directory, cannot install Bayesx in:", bin.dir))
    }
  }

  ## installing on unix systems
  if(ok && unixos) {
    ok <- file.exists(paste(bin.dir, "/bayesxsource.zip", sep = ""))
    if(is.null(source.dir) && !ok) {
      textfun("downloading BayesX sourcecode, may take some time!")
      download.file("http://www.stat.uni-muenchen.de/~bayesx/install/bayesxsource.zip",
        paste(bin.dir, "/bayesxsource.zip", sep = ""))
    }
    if(!is.null(source.dir) && !ok) {
      if(!file.exists(source.dir))
        stop("the source file directory is not existing!")
      if(file.exists(paste(source.dir, "/bayesxsource.zip", sep = ""))) {
        file.copy(paste(source.dir, "/bayesxsource.zip", sep = ""), 
          paste(bin.dir, "/bayesxsource.zip", sep = ""))
      } else {
        options(warn = 0)
        stop(paste("there is no BayesX source file in ", source.dir, "!", sep = ""))
      }
    }
    dir <- paste(bin.dir, "/ ", sep = "")
    ok <- file.exists(paste(bin.dir, "/bayesx", sep = ""))
    if(ok) {
      options(warn = 0)
      stop(paste(bin.dir, "/bayesx ", 
        "does already exist, will not extract bayesxsource.zip!", sep = ""))
    }
    zip.file.extract(file = dir, zipname = "bayesxsource.zip",
      unzip = getOption("unzip"), dir = dir)	
    dir <- paste(strsplit(dir," "," ")[[1]], collapse = "")
    hereisbayesx <- paste(dir, "bayesx/sourcecode", sep = "")
    textfun("start compiling sourcecode, may take some time!")
    owd <- getwd()
    setwd(hereisbayesx)
    ok <- try(system("make", intern = TRUE), silent = TRUE)
    if(length(ok) == 0) {
      options(warn = 0)
      stop("could not compile BayesX sourcecode!")
    }
    ok <- file.exists(base_name)
    if(ok) {
      file.copy(paste(hereisbayesx, "/", base_name, sep = ""), paste(dir, base_name, sep = ""))
      try(system(paste("rm -r -f ", dir, "bayesx", sep = ""), intern = TRUE), silent = TRUE)
      textfun(paste("BayesX was successfully installed in: ", bin.dir, sep = ""))
      if(!is.null(owd))
        try(setwd(owd), silent = TRUE)
      try(Sys.chmod(paste(bin.dir, "/", base_name, sep = ""), mode = "0777"), silent = TRUE)
    } else {
      if(!is.null(owd))
        try(setwd(owd), silent = TRUE)
      options(warn = 0)
      stop("could not install BayesX (maybe see if folder bayesx already exists)!")
    }
  } else {
    ## installing on windows systems
    ok <- file.exists(paste(bin.dir, "/BayesX_windows.exe", sep = ""))
    if(is.null(source.dir) && !ok) {
      textfun("downloading BayesX windows installer, may take some time!")
      download.file("http://www.stat.uni-muenchen.de/~bayesx/install/BayesX_windows.exe",
        paste(bin.dir, "/BayesX_windows.exe", sep = ""), mode = "wb")
    }
    if(!is.null(source.dir) && !ok) {
      if(!file.exists(source.dir))
        stop("the source file directory is not existing!")
      if(file.exists(paste(source.dir, "/BayesX_windows.exe", sep = ""))) {
        file.copy(paste(source.dir, "/BayesX_windows.exe", sep = ""), 
          paste(bin.dir, "/BayesX_windows.exe", sep = ""))
      } else {
        options(warn = 0)
        stop(paste("there is no BayesX windows installer in ", source.dir, "!", sep = ""))
      }
    }
    ok <- file.exists(paste(bin.dir, "/BayesX_windows.exe", sep = ""))
    if(ok) {
      owd <- getwd()
      setwd(bin.dir)
      try(Sys.chmod("BayesX_windows.exe", mode = "0777"), silent = TRUE)
      try(system("BayesX_windows.exe -manual", intern = TRUE), silent = TRUE)
      setwd(owd)
      cat("please supply the path the command line version of BayesX is installed,\n",
        "search folder \'commandline\' in install path of BayesX!\n")
      cmd <- readLines(n = 1)
      ok <- file.exists(cmd, "bayesx.exe")[1L]
      if(ok)
        file.copy(paste(cmd, "/bayesx.exe", sep = ""), paste(bin.dir, "/bayesx.exe", sep = ""))
      else {
        options(warn = 0)
        stop(paste("cannot find the executable bayesx.exe in:", cmd))
      }
      if(ok)
        textfun(paste("BayesX was successfully installed in: ", bin.dir, sep = ""))
    }
  }
  options(warn = 0)

  if(ok)
    return(bin.dir)
  else
    return(paste(bin.dir, "/", base_name, sep = ""))
}


## some fun with text
textfun <- function(text, speed = 0.05)
{
  text <- strsplit(text,""," ")[[1L]]
  Sys.sleep(speed)
  for(w in text) {
    cat(w)
    Sys.sleep(speed)
  }
  Sys.sleep(speed)
  cat("\n")

  return(invisible(NULL))
}


check.install.bayesx <- function(bin = getOption("bayesx.bin"), verbose = FALSE)
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
