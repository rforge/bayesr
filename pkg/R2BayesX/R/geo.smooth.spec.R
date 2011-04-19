geo.smooth.spec <- function(object, dir, prg, data, type)
{
  if(!is.list(object$xt))
    object$xt <- list(object$xt)
  map.name <- help.map.name(deparse(substitute(object, env = .GlobalEnv), 
    backtick = TRUE, width.cutoff = 500L))
  map <- object$xt$map
  if(is.null(map)) {
    if(!is.list(object$xt[[1L]]))
      map <- object$xt
    else {
      map <- NULL
      for(i in 1L:length(object$xt))
        if(inherits(object$xt[[i]], "bnd") || inherits(object$xt[[i]], "list"))
          map <- object$xt[[i]]
    }
    if(is.null(map)) {
      map <- object$xt
      if(is.null(map) || (!is.list(map) && !inherits(map, "bnd")))
        stop("need to supply a bnd file object in argument xt!")
    }
  }
  if(!inherits(map, "bnd"))
    class(map) <- "bnd"
  counter <- NULL
  ok <- TRUE
  if(!missing(dir)) {
    files <- list.files(dir)
    while(ok) {
      mapfile <- paste(map.name, counter, ".bnd", sep = "")
      if(any(grepl(mapfile, files))) {
        if(is.null(counter))
          counter <- 0L
        counter <- counter + 1L
      } else
        ok <- FALSE
    }
    mapfile <- paste(dir, "/", mapfile, sep = "")
    prgfile <- paste(dir, "/", prg, sep = "")
    cat("map", map.name, "\n", file = prgfile, append = TRUE)
    write.bnd(map = map, file = mapfile, replace = TRUE)
    cmd <- paste(map.name, ".infile using ", mapfile, "\n", sep = "")
    cat(cmd, file = prgfile, append = TRUE)
  }	  
  term <- object$term
  if(is.na(object$p.order[1L]))
    object$p.order <- c(3L, 1L)
  if(object$p.order[2L] > 1L) {
    object$p.order[2L] <- 1L
    warning("only random walks of order 1 supported for geosplines, set to default!")
  }
  if(length(object$p.order) < 2L)
    object$p.order <- c(object$p.order, 1L)
  if(object$bs.dim < 0L)
    object$bs.dim <- as.integer(length(map) / 2)
  else {
    if(object$bs.dim >= length(map))
      stop("basis dimension is larger than existing polygons in bnd object!")
  }
  nrknots <- object$bs.dim - object$p.order[1L] + 1L
  if(type == "geokriging") {
    if(!is.null(object$xt$full)) {
      term <- paste(term, "(", type, ",map=", map.name, ",full", sep = "")
      object$xt$full <- NULL
    } else term <- paste(term, "(", type, ",map=", map.name, ",nrknots=", nrknots, sep = "")
  } else term <- paste(term, "(", type, ",map=", map.name, ",nrknots=", nrknots, sep = "")
  term <- paste(do.xt(term, object, "map"), ")", sep = "")
  if(object$by != "NA") {
    if(!is.character(data)) {
      by <- eval(parse(text = object$by), envir = data)
      if(is.factor(by))
        by <- paste(object$by, levels(by), sep = "")
      else
        by <- object$by
    } else by <- object$by
    term <- paste(by, "*", term,sep="")
  }

  return(term)
}
