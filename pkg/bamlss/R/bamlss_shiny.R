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
         uiOutput("select_variables"),
         uiOutput("data_select")
       )
     ),
     ## Show a plot of the generated distribution.
     mainPanel(
       uiOutput("predict"),
       uiOutput("plot_type_button"),
       uiOutput("plot_parameters"),
       uiOutput("plot_button"),
       plotOutput("the_plot", width = '60%'),
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
    if((input$selected_model != "") | input$selected_model != "No 'bamlss' model available!") {
      m <- get(input$selected_model, envir = .GlobalEnv)
      mf <- na.omit(model.frame(m))
      output$select_variables <- renderUI({
        selectInput('variables_selected', 'Select covariates for visualization.',
          names(mf), multiple = TRUE, selectize = TRUE)
      })
      output$data_select <- renderUI({
        if(!is.null(input$variables_selected)) {
          sliders <- lapply(input$variables_selected, function(i) {
            vn <- paste("var_", i, sep = "")
            if(!is.factor(mf[[i]])) {
              x <- pretty(mf[[i]])
              sliderInput(vn, i, min = min(x), max = max(x), value = range(x))
            } else {
              selectInput(vn, i, levels(mf[[i]]))
            }
          })
          rval <- do.call(tagList, c(list(tags$b("Choose covariate values.")), sliders))
        } else rval <- NULL
        rval
      })
    }
  })

  output$predict <- renderUI({
    if(!is.null(input$variables_selected)) {
      rval <- tagList(
        numericInput("ngrid", "Size of grid.", 30, min = 1, max = 1000),
        actionButton("predict_model", label = "Predict."),
        tags$hr()
      )
    } else rval <- NULL
    rval
  })

  observeEvent(input$predict_model, {
    if(!is.null(input$variables_selected)) {
      m <- get(input$selected_model, envir = .GlobalEnv)
      mf <- na.omit(model.frame(m))
      nd <- list()
      for(i in input$variables_selected) {
        vn <- paste("var_", i, sep = "")
        if(is.character(input[[vn]])) {
          nd[[i]] <- factor(input[[vn]], levels = levels(mf[[i]]))
        } else {
          nd[[i]] <- seq(input[[vn]][1], input[[vn]][2], length = as.integer(input$ngrid))
        }
      }
      nd <- expand.grid(nd)
      withProgress(message = "Computing predictions", value = 0, {
        pred <- try(predict(m, newdata = nd, term = input$variables_selected,
          FUN = function(x) { x }), silent = TRUE)
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
    e2 <- input$variables_selected
    e3 <- input$predict_model
    dir("predictions")
  })

  output$plot_type_button <- renderUI({
    rval <- NULL
    if(length(available_predictions())) {
      m <- get(input$selected_model, envir = .GlobalEnv)
      nf <- family(m)$names
      rval <- tagList(
        selectInput('plot_type', 'Select the type of plot.', c("Histogram")),
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
    if(length(available_predictions())) {
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
            hist(pred[[i]], breaks = if(input$nbreaks < 0) "Sturges" else input$nbreaks, freq = FALSE,
              main = i, xlab = "Predictions")
            lines(density(pred[[i]]))
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

