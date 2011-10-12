install.bayesx <-
function(inst.dir = NULL, source.dir = NULL, type = NULL)
{
  warn <- getOption("warn")
  options(warn = -1L)
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

  ## installing/compiling from source
  if(ok && (unixos || (!is.null(type) && type != "binary"))) {
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
    dir <- file.path(bin.dir)
    ok <- file.exists(paste(bin.dir, "/bayesx", sep = ""))
    if(ok) {
      options(warn = 0)
      warning(paste(bin.dir, "/bayesx ", 
        "does already exist, will not extract bayesxsource.zip!", sep = ""))
    }
    #zip.file.extract(file = paste(dir, "/ ", sep = ""), zipname = "bayesxsource.zip",
    #  unzip = getOption("unzip"), dir = paste(dir, "/ ", sep = ""))
    ## begin FIXME
    unzip(file.path(dir, "bayesxsource.zip"), exdir = dir)
    ## end FIXME
    dir <- paste(strsplit(dir," "," ")[[1]], collapse = "")
    hereisbayesx <- file.path(dir, "bayesx/sourcecode")
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
      file.copy(file.path(hereisbayesx, base_name), file.path(dir, base_name))
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
    ## installing the binary on windows systems
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
        warning(paste("cannot find the executable bayesx.exe in:", cmd))
      }
      if(ok)
        textfun(paste("BayesX was successfully installed in: \"", bin.dir, "\"", sep = ""))
    }
  }
  options(warn = warn)

  if(ok)
    return(bin.dir)
  else
    return(paste(bin.dir, "/", base_name, sep = ""))
}

