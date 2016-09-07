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

  dump("bamlss_shiny_ui", file = file.path(dir, "ui.R"), envir = .GlobalEnv)
  dump("bamlss_shiny_server", file = file.path(dir, "server.R"), envir = .GlobalEnv)

  shiny::runApp(dir)
}


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
       selectInput('plot_type', 'Select the type of plot.', c("Histogram")),
       uiOutput("select_parameters"),
       tags$hr(),
       uiOutput("plot_parameters"),
       actionButton("plot_predictions", label = "Plot."),
       plotOutput("the_plot", width = '60%')
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
        numericInput("ngrid", "Size of grid.", 100, min = 1, max = 1000),
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
      pred <- try(predict(m, newdata = nd, term = input$variables_selected), silent = TRUE)
      if(!inherits(pred, "try-error"))
        save(pred, file = "predictions/pred.Rda")
      else warning("Cannot compute prediction, check selected covariates!")
    }
  })

  output$select_parameters <- renderUI({
    rval <- NULL
    if(!is.null(input$variables_selected)) {
      m <- get(input$selected_model, envir = .GlobalEnv)
      nf <- family(m)$names
      if(length(nf) > 1) {
        rval <- selectInput('selected_parmeters', 'Select distribution parameter.',
          nf, selected = nf[1], multiple = TRUE, selectize = TRUE)
      }
    }
    rval
  })

  output$plot_parameters <- renderUI({
    if(input$plot_type == "Histogram") {
      rval <- tagList(numericInput("nbreaks", "Number of breaks.", -1, min = 1, max = 1000, step = 1))
    } else rval <- NULL
    rval
  })

  observeEvent(input$plot_predictions, {
    if(file.exists("predictions/pred.Rda")) {
      output$the_plot <- renderPlot({
        load("predictions/pred.Rda")
        par(mfrow = n2mfrow(length(pred)), mar = c(4.1, 4.1, 4.1, 0.1))
        nf <- names(pred)
        if(length(nf) > 1)
          nf <- input$selected_parmeters
        for(i in nf) {
          hist(pred[[i]], breaks = if(input$nbreaks < 0) "Sturges" else input$nbreaks, freq = FALSE,
            main = i, xlab = "Predictions")
        }
      })
    }
  })
}

