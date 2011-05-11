bayesx_prgfile <-
function(x, model = 1L)
{
  bayesx.prg <- NULL
  if(inherits(x, "bayesx")) {
    bayesx.prg <- x[[model]]$bayesx.prg$prg
      if(!is.null(bayesx.prg))
        cat(bayesx.prg)
  }

  return(invisible(bayesx.prg))
}

