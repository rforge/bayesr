read.bayesx.model.output <-
function(dir, model.name)
{
  if(is.null(dir))
    stop("no directory specified!")
  if(is.null(model.name))
    stop("no model name specified!")
  files <- dir.files <- list.files(dir)
  if(!any(grep(model.name,files)))
    stop(paste("no model results existing for ", model.name, "!", sep = ""))
  else {
    rval <- list()
    files <- grep(model.name, files, value = TRUE, fixed = TRUE)
    filep <- grep(paste(model.name, ".", sep = ""), files, value = TRUE, fixed = TRUE)  
    files <- c(grep(paste(model.name, "_", sep = ""), files, value = TRUE, fixed = TRUE), filep)  

    ## search for data
    data <- NULL
    if(length(i <- grep(paste(model.name, ".data.raw", sep = ""), files)))
      data <- as.matrix(read.table(paste(dir, "/", files[i], sep = ""), header = TRUE))
    for(char in c("_predict.raw", "_predictmean.raw", "_predict.res"))
      if(length(i <- grep(char, files)))
        data <- cbind(data, as.matrix(read.table(paste(dir, "/", files[i], sep = ""), header = TRUE)))

    ## set response and predictor
    N <- NA
    response <- eta <- residuals <- NULL
    if(!is.null(data)) {
      data[data == "."] <- "NA"
      mode(data) <- "numeric"
      data <- as.data.frame(data)
      data$intnr <- NULL
      dn <- unique(names(data))
      data <- data[dn]
      if(is.null(response))
        response <- data[[1L]]
      if(is.factor(response)) 
        response <- as.numeric(as.character(response))
      N <- length(response)
      rval$fitted.values <- eta <- get.eta(data)
      rval$residuals <- response - eta
      rval$response <- response
    }
    ## get smooth and random effects
    rval <- c(rval, find.smooth.random(dir, files, data, response, eta, model.name))

    ## get fixed effects
    rval <- find.fixed.effects(dir, files, data, response, eta, model.name, rval)

    ## get scale estimate
    rval$variance <- get.scale(files, dir)

    ## search for other results
    model.results <- mf <- mf2 <- NULL
    method <- ""
    if(any(grep("deviance.raw", files))) {
      mf <- grep("deviance.raw", files, value = TRUE)
      mf <- read.table(paste(dir, "/", mf, sep = ""), header = TRUE)
      pd <- mf$unstandardized_deviance[length(mf$unstandardized_deviance) - 1L]
      DIC <- mf$unstandardized_deviance[length(mf$unstandardized_deviance)]
      mf <- list(DIC = DIC, pd = pd)
    }
    if(any(grep(".tex",files))) {
      sm <- readLines(paste(dir, "/", grep(".tex", files, value = TRUE), sep = ""))
      model.results <- search.bayesx.tex(sm)
      method <- model.results$method
    }
    if(any(grep("modelfit.raw",files))) {
      mf <- grep("modelfit.raw", files, value = TRUE)
      mf <- chacol(read.table(paste(dir, "/", mf, sep = ""), header = TRUE))
      mf <- mf2 <- as.list(mf)
    }
    if(!is.null(model.results)) {
      if(!is.null(mf2)) {
        n1 <- names(mf2)
        n2 <- names(model.results)
        for(i in 1L:length(mf2))
          if(!n1[i] %in% n2)
            eval(parse(text = paste("model.results$", n1[i], "<- mf2[[i]]", sep = "")))
      }
      mf <- model.results
    } else mf <- c(list(method = method, N = N), mf)
    if(mf$method == " ")
      mf$method <- "NA"
    if(!is.null(mf$logLik))
      mf$logLik <- mf$logLik/(-2)
    rval$model.fit <- mf

    ## reordering and naming
    if((info <- paste(model.name, ".terms.info", sep = "")) %in% dir.files)
      rval$effects <- term.reorder(rval$effects, paste(dir, "/", info, sep = ""))

    ## search for additional info
    rval$model.fit <- smi(info, rval$model.fit)

    ## reformate output
    rval <- bayesx.reformate(rval)

    return(rval)
  }
}

