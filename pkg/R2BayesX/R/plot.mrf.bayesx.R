plot.mrf.bayesx <-
function(x, diagnostics = FALSE, ...)
{
  if(!is.null(x)) {
		args <- list(...)
    if(diagnostics == FALSE) {
      if(inherits(x, "random.bayesx")) {
        cx <- colnames(x)
        if(!is.null(args$total) && args$total && any(grep("_tot", cx))) {
          xattr <- attributes(x)
          id <- c(1L, grep("_tot", cx))
          x <- x[,id]
          for(na in names(xattr))
            if(na != "dim" && na != "dimnames" && na != "names")
              attr(x, na) <- eval(parse(text = paste("xattr$", na, sep = "")))
          args$total <- NULL
        } else {
          xattr <- attributes(x)
          id <- 1L:ncol(x)
          id <- id[!grepl("tot_sim", cx) & !grepl("_tot", cx) & !grepl("_sim", cx)]
          x <- x[,id]
          for(na in names(xattr))
            if(na != "dim" && na != "dimnames" && na != "names")
              attr(x, na) <- eval(parse(text = paste("xattr$", na, sep = "")))
        }
      }
      if(!is.null(attr(x, "map.name"))) {
        args$map <- try(eval(parse(text = attr(x, "map.name")), envir = globalenv()), silent = TRUE)
        if(class(args$map) == "try-error")
          args$map <- NULL
      }
      if(is.null(args$map)) {
        args$x <- x
        do.call("plotblock", args)
      } else {
        args$x <- x
        do.call("plotmap", args)
      }
    } else coeffun(x, args, diagnostics)
  } else warning("there is nothing to plot!")

  return(invisible(NULL))
}

