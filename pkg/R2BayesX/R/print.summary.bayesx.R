print.summary.bayesx <- function(x, digits = max(3, getOption("digits") - 3),
  signif.stars = getOption("show.signif.stars"), ...)
{
  n <- length(x)
  ncheck <- n > 1L
  digits <- attr(x, "digits")
  for(i in 1L:n) {
    if(ncheck)
      cat("###", i, "\n")
    .print_summary_bayesx(x[[i]], digits = digits, signif.stars = signif.stars, ...)
  }
  if(ncheck) {
    cat("###\n")
    cat("Object contains of", n, "models\n")
  }

  return(invisible(NULL))
}


.print_summary_bayesx <- function(x, digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"), ...)
{
  if(!is.null(x$model.fit))
    if(!is.null(x$model.fit$model.name))
      if(length(grep("_hlevel", x$model.fit$model.name))) {
        hlevel <- splitme(strsplit(x$model.fit$model.name, "_hlevel")[[1]][2])
        go <- TRUE
        hl <- NULL
        for(i in 1:length(hlevel)) {
          if(hlevel[i] == "_")
            go <- FALSE
          if(go)
            hl <- c(hl, hlevel[i])
        }
        hlevel <- as.integer(resplit(hl))
        if(hlevel > 1)
          cat("Hierarchical random effects model results: stage", hlevel, "\n")
        else {
          cat("Main effects model results: stage", hlevel, "\n")
          cat("---\n")
        }
      }
  if(!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
  } else {
    if(!is.null(x$model.fit$formula)) {
      cat("Formula:\n")
      print(as.formula(x$model.fit$formula))
    }
  }
  liner <- "---"
  fc <- FALSE
  if(!is.null(x$fixed.effects)) {
    fc <- TRUE
    if(nrow(x$fixed.effects) < 2L) {
      if(all(x$fixed.effects[1L,] == 0))
        fc <- FALSE
    } else {
      if(all(x$fixed.effects[1L,] == 0)) {
        m <- ncol(x$fixed.effects)
        nc <- colnames(x$fixed.effects)
        nr <- rownames(x$fixed.effects)[2L:nrow(x$fixed.effects)]
        x$fixed.effects <- matrix(x$fixed.effects[2L:nrow(x$fixed.effects),], ncol = m)
        colnames(x$fixed.effects) <- nc
        rownames(x$fixed.effects) <- nr
      }
    }
    x$fixed.effects <- round(x$fixed.effects, digits)
  }
  if(fc || (!is.null(x$smooth.hyp))) {
    cat(liner, "\n")
    cat("Fixed effects estimation results:\n")
    cat("-\n")
  }
  if(fc) {
    cat("Parametric Coefficients:\n")
    printCoefmat(x$fixed.effects)
  }
  if(!is.null(x$smooth.hyp)) {
    if(fc)
      cat("-\n")
    cat("Smooth terms variances:\n")
    ls <- ncol(x$smooth.hyp)
    terms <- colnames(x$smooth.hyp)	
    rn <- rownames(x$smooth.hyp)
    x$smooth.hyp <- round(x$smooth.hyp, digits)
    printCoefmat(x$smooth.hyp)
  }
  cat(liner, "\n")
  if(!is.null(x$random.hyp)) {
    cat("Random effects variances:\n")
    x$random.hyp <- round(x$random.hyp, digits)
    printCoefmat(x$random.hyp)		
    cat(liner, "\n")
  }		
  if(!is.null(x$model.fit)) {
    if(x$model.fit$method == "MCMC") {
      if(!is.null(x$variance)) {
        cat("Scale estimate:\n")
        x$variance <- round(x$variance, digits)
        printCoefmat(x$variance)
        cat(liner, "\n")
      }
    } else {
      if(!is.null(x$variance)) {
        cat("Scale estimate:", round(as.numeric(x$variance)[1], digits), "\n")
        cat(liner, "\n")
      }
    }
    mfn <- names(x$model.fit)
    step <- 5L
    for(i in 1L:length(mfn)) {
      if(!mfn[i] %in% c("model.name", "formula", "step.final.model")) { 
        if(!is.null(x$model.fit[[i]]) && x$model.fit[[i]] != "") {
          if(i < step)
            cat(mfn[i], "=", x$model.fit[[i]], " ")
          if(i == step) {
            if(i != length(mfn))
              cat("\n")
            cat(mfn[i], "=", x$model.fit[[i]], " ")
            step <- step + step
          }
        }
      }
    }
    cat("\n")
  }
  if(!is.null(x$model.fit$step.final.model)) {
    cat(liner,"\n")
    cat("Stepwise final model:\n")
    cat("-\n")
    cat(x$model.fit$step.final.model, "\n")
  }

  return(invisible(NULL))
}
