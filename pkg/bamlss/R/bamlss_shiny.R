bamlss_shiny <- function(dir = NULL)
{
  stopifnot(requireNamespace("shiny"))
  stopifnot(requireNamespace("tools"))

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
         uiOutput("select_terms"),
         uiOutput("data_select")
       )
     ),
     ## Show a plot of the generated distribution.
     mainPanel(
       uiOutput("predict"),
       uiOutput("plot_type_button"),
       uiOutput("plot_parameters"),
       uiOutput("plot_button"),
       plotOutput("the_plot", height = "800px", width = "600px"),
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
        bamlss:::all.labels.formula(x$formula)
      })
      tn <- unique(do.call("c", tn))
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
              return(selectInput(vn, i, levels(mf[[i]])))
            }
          })
          rval <- do.call(tagList, c(list(tags$b("Choose data range.")), sliders))
        } else rval <- NULL
        rval
      })
    }
  })

  output$predict <- renderUI({
    if(!is.null(input$terms_selected)) {
      rval <- tagList(
        numericInput("ngrid", "Size of grid.", 30, min = 1, max = 1000),
        selectInput('predict_fun', 'Select function applied on samples.',
          c("None", "Mean", "Median", "Sd", "Var", "95% CI", "99% CI", "c95"), selected = "c95"),
        actionButton("predict_model", label = "Predict."),
        tags$hr()
      )
    } else rval <- NULL
    rval
  })

  observeEvent(input$predict_model, {
    if(!is.null(input$terms_selected)) {
      m <- get(input$selected_model, envir = .GlobalEnv)
      mf <- na.omit(model.frame(m))
      vars <- NULL
      for(j in input$terms_selected)
        vars <- c(vars, all.vars(as.formula(paste("~", j))))
      vars <- unique(vars)
      for(j in seq_along(vars)) {
        if(!any(vars[j] %in% names(mf)))
          vars[j] <- grep(vars[j], names(mf), fixed = TRUE, value = TRUE)
      }
      nd <- list()
      for(i in vars) {
        vn <- paste("var_", i, sep = "")
        if(is.character(input[[vn]])) {
          nd[[i]] <- factor(input[[vn]], levels = levels(mf[[i]]))
        } else {
          nd[[i]] <- seq(input[[vn]][1], input[[vn]][2], length = as.integer(input$ngrid))
        }
      }
      nd <- expand.grid(nd)
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
      withProgress(message = "Computing predictions", value = 0, {
        pred <- try(predict(m, newdata = nd, term = input$terms_selected, FUN = FUN), silent = TRUE)
      })
      if(file.exists("predictions/pred.Rda"))
        file.remove("predictions/pred.Rda")
      if(!inherits(pred, "try-error")) {
        save(pred, file = "predictions/pred.Rda")
      } else {
        warning("Cannot compute prediction, check selected covariates!")
      }
    }
  })

  available_predictions <- reactive({
    e1 <- input$load_model
    e2 <- input$terms_selected
    e3 <- input$predict_model
    dir("predictions")
  })

  output$plot_type_button <- renderUI({
    rval <- NULL
    if(length(available_predictions())) {
      m <- get(input$selected_model, envir = .GlobalEnv)
      nf <- family(m)$names
      rval <- tagList(
        selectInput('plot_type', 'Select the type of plot.', c("Histogram", "Effect"), selected = "Histogram"),
        if(length(nf) > 1) {
          selectInput('selected_parameters', 'Select distribution parameter.',
            nf, selected = nf[1], multiple = TRUE, selectize = TRUE)
        } else NULL,
        tags$hr()
      )
    }
    rval
  })

  output$plot_parameters <- renderUI({
    rval <- NULL
    if(length(available_predictions()) & length(input$plot_type)) {
      if(input$plot_type == "Histogram") {
        rval <- tagList(numericInput("nbreaks", "Number of breaks.", -1, min = 1, max = 1000, step = 1))
      }
    }
    rval
  })

  output$plot_button <- renderUI({
    if(length(available_predictions())) {
      rval <- actionButton("plot_predictions", label = "Plot.")
    } else rval <- NULL
    rval
  })

  observeEvent(input$plot_predictions, {
    if(length(available_predictions())) {
      output$the_plot <- renderPlot({
        if(length(available_predictions()) & file.exists("predictions/pred.Rda")) {
          load("predictions/pred.Rda")
          par(mfrow = n2mfrow(length(pred)), mar = c(4.1, 4.1, 4.1, 0.1))
          nf <- names(pred)
          if(length(nf) > 1)
            nf <- input$selected_parameters
          for(i in nf) {
            if(input$plot_type == "Histogram") {
              rdens <- density(pred[[i]])
              rh <- hist(pred[[i]], plot = FALSE)
              args <- list()
              args$ylim <- c(0, max(c(rh$density, rdens$y)))
              args$freq <- FALSE
              args$x <- pred[[i]]
              args$ylab <- "Density"
              args$xlab <- "Predictions"
              args$main <- i
              args$breaks <- if(input$nbreaks < 0) "Sturges" else input$nbreaks
              ok <- try(do.call("hist", args))
              if(!inherits(ok, "try-error"))
                lines(rdens)
              box()
            }
            if(input$plot_type == "Effect") {
              m <- get(input$selected_model, envir = .GlobalEnv)
              mf <- na.omit(model.frame(m))
              vars <- NULL
              for(j in input$terms_selected)
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
                  nd[[j]] <- seq(input[[vn]][1], input[[vn]][2], length = as.integer(input$ngrid))
                }
              }
              nd <- expand.grid(nd)
              if(length(vars) < 2)
                plot2d(pred[[i]] ~ nd[[vars]], xlab = vars, ylab = "Predictions")
              if(length(vars) < 3)
                plot3d(pred[[i]] ~ nd[[vars[1]]] + nd[[vars[2]]], xlab = vars[1], ylab = vars[2], zlab = "Predictions")
            }
          }
        }
      })
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

