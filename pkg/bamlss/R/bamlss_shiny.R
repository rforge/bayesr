bamlss_shiny <- function(dir = NULL)
{
  owd <- getwd()
  if(is.null(dir)) {
    dir.create(dir <- tempfile())
    on.exit(unlink(dir))
  }
  dir <- path.expand(dir)
  if(!file.exists(file.path(dir, "www")))
    dir.create(file.path(dir, "www"))
  if(!file.exists(file.path(dir, "predictions")))
    dir.create(file.path(dir, "predictions"))

  ## ui.R and server.R
bamlss_shiny_ui <- function(...) {
  fluidPage(
    titlePanel("Model Visualizer"),
    sidebarLayout(
     sidebarPanel(
       tabPanel("Choose a model",
         br(),
         uiOutput("select_model"),
         actionButton("load_model", label = "Load model."),
         tags$hr(),
         uiOutput("select_intercept"),
         uiOutput("select_terms"),
         uiOutput("data_select")
       )
     ),
     ## Show a plot of the generated distribution.
     mainPanel(
       uiOutput("predict_parameters"),
       uiOutput("set_plot"),
       plotOutput("make_plot"),
       tags$br(),
       tags$br()
     )
   ))
}


bamlss_shiny_server <- function(input, output, session)
{
  available_models <- reactive({
    models <- NULL
    for(m in ls(envir = .GlobalEnv)) {
      tm <- get(m, envir = .GlobalEnv)
      if(inherits(tm, "bamlss"))
        models <- c(models, m)
    }
    if(!is.null(input$selected_model)) {
      if(input$selected_model != "") {
        if(input$selected_model %in% models) {
          i <- which(models == input$selected_model)
          models <- c(models[i], models[-i])
        }
      }
    }
    if(is.null(models))
      models <- "No 'bamlss' model available!"
    return(models)
  })

  output$select_model <- renderUI({
    selectInput('selected_model', 'Select model for visualization.',
      available_models())
  })

  observeEvent(input$load_model, {
    if(file.exists("predictions/pred.Rda"))
      file.remove("predictions/pred.Rda")
    if((input$selected_model != "") & !grepl("model available", input$selected_model)) {
      m <- get(input$selected_model, envir = .GlobalEnv)
      mf <- model.frame(m)
      tn <- lapply(formula(m), function(x) {
        all.labels.formula(x$formula)
      })
      tn <- unique(do.call("c", tn))
      output$select_intercept <- renderUI({
        selectInput("intercept", "Include intercept?", c("yes", "no"), selected = "no")
      })
      output$select_terms <- renderUI({
        selectInput('terms_selected', 'Select terms for visualization.',
          tn, multiple = TRUE, selectize = TRUE)
      })
      output$data_select <- renderUI({
        if(!is.null(input$terms_selected)) {
          vars <- NULL
          for(j in input$terms_selected)
            vars <- c(vars, all.vars(as.formula(paste("~", j))))
          vars <- unique(vars)
          for(j in seq_along(vars)) {
            if(!any(vars[j] %in% names(mf)))
              vars[j] <- grep(vars[j], names(mf), fixed = TRUE, value = TRUE)
          }
          sliders <- lapply(vars, function(i) {
            vn <- paste("var_", i, sep = "")
            if(!is.factor(mf[[i]])) {
              x <- pretty(mf[[i]])
              return(sliderInput(vn, i, min = min(x), max = max(x), value = range(x)))
            } else {
              return(selectInput(vn, i, levels(mf[[i]]), multiple = TRUE, selectize = TRUE, selected = levels(mf[[i]])[1]))
            }
          })
          rval <- do.call(tagList, c(list(tags$b("Choose data range/levels.")), sliders))
        } else rval <- NULL
        rval
      })
    }
  })

  output$predict_parameters <- renderUI({
    if(!is.null(input$selected_model)) {
      m <- get(input$selected_model, envir = .GlobalEnv)
      nf <- c(family(m)$names, family(m)$moment)
      nsamps <- nrow(m$samples[[1]]) - 1L
      rval <- tagList(
        numericInput("ngrid", "Size of grid.", 30, min = 1, max = 1000),
        numericInput("nsamps", "Number of samples.", nsamps, min = 1, max = nsamps),
        selectInput('predict_fun', 'Select function applied on samples.',
          c("None", "Mean", "Median", "Sd", "Var", "95% CI", "99% CI", "c95"), selected = "c95"),
        if(length(nf) > 1) {
          selectInput('selected_parameters', 'Select distribution parameter/moment.',
            nf, selected = nf[1])
        } else NULL,
        selectInput('plot_type', 'Select the type of plot.',
          c("Histogram", "Histogram 2d", "Effect", "Image"), selected = "Histogram"),
        conditionalPanel(condition = "input.plot_type == 'Histogram'",
          numericInput("nbreaks", "Number of breaks.", -1, min = 1, max = 1000, step = 1)),
        tags$hr()
      )
    } else rval <- NULL
    rval
  })

  model_terms <- reactive({
    e1 <- input$intercept
    input$terms_selected
  })

  newdata <- reactive({
    nd <- NULL
    if(!is.null(input$intercept))
      nd <- input$intercept
    if(!is.null(model_terms()) & !is.null(input$ngrid)) {
      m <- get(input$selected_model, envir = .GlobalEnv)
      mf <- na.omit(model.frame(m))
      vars <- NULL
      for(j in model_terms())
        vars <- c(vars, all.vars(as.formula(paste("~", j))))
      vars <- unique(vars)
      for(j in seq_along(vars)) {
        if(!any(vars[j] %in% names(mf)))
          vars[j] <- grep(vars[j], names(mf), fixed = TRUE, value = TRUE)
      }
      nd <- list()
      all_NULL <- FALSE
      for(i in vars) {
        vn <- paste("var_", i, sep = "")
        if(!is.null(input[[vn]])) {
          if(is.character(input[[vn]])) {
            nd[[i]] <- factor(input[[vn]], levels = levels(mf[[i]]))
          } else {
            nd[[i]] <- if(identical(input[[vn]][1], input[[vn]][2])) {
              input[[vn]][1]
            } else {
              seq(input[[vn]][1], input[[vn]][2], length = as.integer(input$ngrid))
            }
          }
          all_NULL <- c(all_NULL, FALSE)
        } else {
          all_NULL <- c(all_NULL, TRUE)
        }
      }
      nd <- if(!all(all_NULL)) expand.grid(nd) else NULL
    }
    nd
  })

  predict_model <- reactive({
    if(!is.null(nd <- newdata()) & !is.null(input$predict_fun)) {
      m <- get(input$selected_model, envir = .GlobalEnv)
      FUN <- switch(input$predict_fun,
        "None" = function(x) { x },
        "Mean" = function(x) { mean(x, na.rm = TRUE) },
        "Median" = function(x) { median(x, na.rm = TRUE) },
        "Sd" = function(x) { sd(x, na.rm = TRUE) },
        "Var" = function(x) { var(x, na.rm = TRUE) },
        "c95" = c95,
        "95% CI" = function(x) { quantile(x, probs = c(0.025, 0.975), na.rm = TRUE) },
        "99% CI" = function(x) { quantile(x, probs = c(0.005, 0.995), na.rm = TRUE) }
      )
      if(is.character(nd)) {
        nd <- NULL
      } else {
        if(nrow(nd) < 1)
          nd <- NULL
      }
      if(is.null(model_terms()) & input$intercept == "no") {
        pred <- NULL
      } else {
        withProgress(message = "Generating visualization data ...", value = 0, {
          mt <- model_terms()
          if(is.null(mt) & (input$intercept == "yes"))
            mt <- "Intercept"
          pred <- try(predict(m, newdata = nd, term = mt,
            FUN = FUN, intercept = input$intercept == "yes", nsamps = input$nsamps), silent = TRUE)
        })
        if(inherits(pred, "try-error"))
          pred <- NULL
      }
    } else pred <- NULL
    pred
  })

  output$make_plot <- renderPlot({
    if(!is.null(pred <- predict_model()) & length(input$selected_parameters) & length(input$ngrid)) {
      nf <- names(pred)
      if(length(nf) > 1)
        nf <- input$selected_parameters
      if(input$plot_type != "Image")
        par(mfrow = n2mfrow(length(nf)))
      is_hist <- input$plot_type %in% c("Histogram", "Histogram 2d")
      for(i in nf) {
        if(input$plot_type == "Histogram") {
          par(mar = c(4.1, 4.1, if(length(nf) > 1) 4.1 else 1.1, 0.1))
          rdens <- density(pred[[i]])
          rh <- hist(pred[[i]], plot = FALSE)
          args <- list()
          args$ylim <- c(0, max(c(rh$density, rdens$y)))
          args$freq <- FALSE
          args$x <- pred[[i]]
          args$ylab <- "Density"
          args$xlab <- "Predictions"
          args$main <- if(length(nf) > 1) i else NA
          args$col <- "lightgray"
          args$breaks <- if(input$nbreaks < 0) "Sturges" else input$nbreaks
          ok <- try(do.call("hist", args))
          if(!inherits(ok, "try-error"))
            lines(rdens)
          box()
        }
        if(input$plot_type != "Histogram" & !is.null(model_terms())) {
          m <- get(input$selected_model, envir = .GlobalEnv)
          mf <- na.omit(model.frame(m))
          vars <- NULL
          for(j in model_terms())
            vars <- c(vars, all.vars(as.formula(paste("~", j))))
          vars <- unique(vars)
          for(j in seq_along(vars)) {
            if(!any(vars[j] %in% names(mf)))
              vars[j] <- grep(vars[j], names(mf), fixed = TRUE, value = TRUE)
          }
          nd <- list()
          for(j in vars) {
            vn <- paste("var_", j, sep = "")
            if(is.character(input[[vn]])) {
              nd[[j]] <- factor(input[[vn]], levels = levels(mf[[j]]))
            } else {
              nd[[j]] <- if(identical(input[[vn]][1], input[[vn]][2])) {
                input[[vn]][1]
              } else {
                seq(input[[vn]][1], input[[vn]][2], length = as.integer(input$ngrid))
              }
            }
          }
          nd <- expand.grid(nd)
          dim <- apply(nd, 2, function(x) {
            all(x == x[1])     
          })
          vars <- vars[!dim]
          if(length(vars) == 1 & (input$plot_type == "Effect")) {
            par(mar = c(4.1, 4.1, if(length(nf) > 1) 4.1 else 1.1, 0.1))
            plot2d(pred[[i]] ~ nd[[vars]], xlab = vars, ylab = "Predictions", ylim = range(pred[[i]]),
              col.lines = if(input$predict_fun == "None") rgb(0.1, 0.1, 0.1, alpha = 0.1) else "black",
              main = if(length(nf) > 1) i else NULL,
              fill.select = if(input$predict_fun == "c95") c(0, 1, 0, 1) else NULL,
              scheme = if(input$predict_fun == "c95") 2 else 1, grid = 50)
          }
          if(length(vars) == 2) {
            if(input$plot_type %in% c("Effect", "Image")) {
              par(mar = if(input$plot_type == "Image") {
                  c(4.1, 4.1, if(length(nf) > 1) 4.1 else 1.1, 0.1)
                } else c(1.1, 0.1, if(length(nf) > 1) 4.1 else 0.1, 0.1))
              c.select <- if(is.null(dim(pred[[i]]))) 1 else ncol(pred[[i]])
              plot3d(pred[[i]] ~ nd[[vars[1]]] + nd[[vars[2]]], xlab = vars[1], ylab = vars[2],
                zlab = "Predictions", c.select = 1:c.select,
                border = if(input$predict_fun == "None") rgb(0.1, 0.1, 0.1, alpha = 0.01) else NULL,
                image = input$plot_type == "Image", main = if(length(nf) > 1) i else NULL,
                legend = length(nf) < 2, grid = input$ngrid, ticktype = "detailed",
                contour = input$plot_type == "Image")
            }
            if(input$plot_type %in% c("Histogram 2d")) {
              
            }
          }
        }
      }
    }
  })
}
  
  nenv <- new.env()
  on.exit(unlink(nenv), add = TRUE)

  assign("bamlss_shiny_ui", bamlss_shiny_ui, envir = nenv)
  assign("bamlss_shiny_server", bamlss_shiny_server, envir = nenv)

  dump("bamlss_shiny_ui", file = file.path(dir, "ui.R"), envir = nenv)
  dump("bamlss_shiny_server", file = file.path(dir, "server.R"), envir = nenv)

  shiny::runApp(dir)
}

