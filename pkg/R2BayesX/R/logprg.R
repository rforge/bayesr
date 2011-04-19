bayesx_prgfile <- function(x, model = 1L)
{
  bayesx.prg <- NULL
  if(inherits(x, "bayesx")) {
    bayesx.prg <- x[[model]]$bayesx.prg$prg
      if(!is.null(bayesx.prg))
        cat(bayesx.prg)
  }

  return(invisible(bayesx.prg))
}


bayesx_logfile <- function(x, model = 1L)
{
  bayesx.log <- NULL
  if(inherits(x, "bayesx")) {
    bayesx.log <- x[[model]]$bayesx.run$log
      if(!is.null(bayesx.log)) {
        if(!is.character(bayesx.log))
          print(bayesx.log)
        else
          writeLines(bayesx.log)
      }
  }

  return(invisible(bayesx.log))
}
