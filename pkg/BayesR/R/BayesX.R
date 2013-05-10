#######################################
## (1) BayesX model fitting wrapper. ##
#######################################
bayesx2 <- function(formula, family = gaussian(), data = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.fail, contrasts = NULL,
  cores = NULL, combine = TRUE, ...)
{
  pm <- match.call(expand.dots = TRUE)
  pm$parse.input <- as.name("parse.input.bayesr")
  pm$transform <- as.name("tranformBayesX")
  pm$setup <- as.name("setupBayesX")
  pm$sampler <- as.name("samplerBayesX")
  pm$results <- as.name("resultsBayesX")
  pm[[1]] <- as.name("bayesr")
  pm <- eval(pm, parent.frame())
  pm
}


####################################
## (2) BayesX specific functions. ##
####################################
tranformBayesX <- function(x, ...)
{
  if(inherits(x, "list") & !("smooth" %in% names(x))) {
    for(j in seq_along(x))
      x[[j]] <- tranformBayesX(x[[j]])
  } else x <- randomize(x)
  x
}


setupBayesX <- function(x, ...)
{
  args <- list(...)
  if(is.null(dir <- args$dir))
    dir.create(dir <- tempfile())
  name <- if(is.null(args$name)) "bayesr" else args$name
  data.name <- if(is.null(args$data.name)) "data" else args$data.name

  if(inherits(x, "list") & !("smooth" %in% names(x))) {
    n <- length(x)
  } else {
    x <- list(x)
    n <- 1
  }

  prg <- paste('% usefile', file.path(dir, paste(name, 'prg', sep = '.')))
  prg <- c(prg, paste('logopen using', file.path(dir, paste(name, 'log', sep = '.'))))
  prg <- c(prg, paste('mcmcreg', name))

#dataset d1  
#d1.infile using /tmp/Rtmpybbo5u/bayesx/bayesx.estim_hlevel1_MAIN_REGRESSION.data.raw 
# 
#b.outfile = /tmp/Rtmpybbo5u/bayesx/bayesx.estim_hlevel1_MAIN_REGRESSION 
 
  eqs <- list()
  for(j in 1:n) {
    if("response" %in% names(x[[j]])) {
      prg2 <- paste("b.hregress", x[[j]]$response)
      if(all(x[[j]]$X[, 1] == 1))
        prg2 <- paste(prg2, "=", "const")
      if(length(x[[j]]$pterms))
        prg2 <- paste(prg2, "+", paste(x[[j]]$pterms, collapse = " + "))
    } else {

    }
    
  }
print(prg2)
stop()
}


samplerBayesX <- function(x, ...) { x }
resultsBayesX <- function(x, ...) { x }
