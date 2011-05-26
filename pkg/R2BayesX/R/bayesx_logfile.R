bayesx_logfile <-
function(x, model = 1L)
{
  x <- get.model(x, model)
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

