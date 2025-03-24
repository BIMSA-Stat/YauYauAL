ui <- fluidPage(
  theme = shinytheme("flatly"),
  tags$head(
    tags$style(HTML("
      /* Overall page padding */
      body { padding-top: 20px; }
      /* Sidebar styling */
      .well {
        background-color: #f7f7f7;
        border: 1px solid #ddd;
        border-radius: 4px;
        padding: 15px;
      }
      /* Button styling */
      .btn-primary {
          background-color: #337ab7;
          border-color: #2e6da4;
          font-size: 16px;
          padding: 8px 20px;
      }
      .btn-success {
          background-color: #5cb85c;
          border-color: #4cae4c;
          font-size: 16px;
          padding: 8px 20px;
      }
      /* Title styling */
      .shiny-title-panel { margin-bottom: 20px; }
      .nav-tabs > li > a {
          font-size: 16px;
          font-weight: bold;
      }
    "))
  ),
  titlePanel("YauYauAL 0.1.0"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      # Arrange controls in fluid rows and columns
      fluidRow(
        column(6,
               numericInput("Dim", "Dimensions:", value = 2, min = 1)
        ),
        column(6,
               numericInput("T", "T:", value = 20, min = 0.1)
        )
      ),
      fluidRow(
        column(6,
               numericInput("Dt", "Dt:", value = 0.001, min = 0.0001)
        ),
        column(6,
               numericInput("Ds", "Ds:", value = 0.5, min = 0.01)
        )
      ),
      textAreaInput("f_expr", "State f(x):",
                    value = "function(x) { return(c(cos(x[1]), cos(x[2]))) }",
                    rows = 2),
      textAreaInput("h_expr", "Observation h(x):",
                    value = "function(x) { return(c(x[1]^3, x[2]^3)) }",
                    rows = 2),
      fluidRow(
        column(6,
               numericInput("seed", "Seed:", value = 42, min = 1, step = 1)
        )
      ),
      fluidRow(
        column(12,
               actionButton("simulate", "Run Simulation", class = "btn-primary")
        )
      ),
      br(),
      fluidRow(
        column(12,
               actionButton("estimate", "Estimate Iu", class = "btn-success")
        )
      ),
      br(),
      fluidRow(
        column(12,
               downloadButton("downloadPlot", "Save Plot as PDF")
        )
      )
    ),
    mainPanel(
      width = 9,
      # Main tabset with updated tab names
      tabsetPanel(
        id = "mainTab",
        tabPanel("State Trajectories", uiOutput("state_plots")),
        tabPanel("Observations", uiOutput("obs_plots"))
      )
    )
  )
)

server <- function(input, output, session) {

  # Reactive values to store simulation and Iu estimation results
  simData <- reactiveValues(x = NULL, y = NULL, Iu = NULL)

  # Run simulation when "Run Simulation" button is clicked
  observeEvent(input$simulate, {
    Dim   <- input$Dim
    T_val <- input$T
    Dt    <- input$Dt
    Ds    <- input$Ds

    # Set Dtau = 5 * Dt
    Dtau    <- 5 * Dt
    Ntau    <- as.integer(T_val / Dtau)
    NtNtau  <- as.integer(T_val / Dt)

    # Parse the functions from the text inputs
    f <- eval(parse(text = input$f_expr))
    h <- eval(parse(text = input$h_expr))
    df <- generate_derivative(f)

    # Simulate state and observation data
    simResult <- Simulate_State_Obser(Dt, Ntau, NtNtau, f, h, Dim, seed = input$seed)
    simData$x <- simResult$x
    simData$y <- simResult$y
    simData$Iu <- NULL  # Clear previous Iu result

    cat("Simulation completed.\n")
  })

  # Estimate Iu with progress indicator when "Estimate Iu" is clicked
  observeEvent(input$estimate, {
    req(simData$x, simData$y)
    withProgress(message = "Computing Iu...", detail = "Preparing parameters", value = 0, {
      Dim   <- input$Dim
      T_val <- input$T
      Dt    <- input$Dt
      Ds    <- input$Ds
      Dtau  <- 5 * Dt
      Nt    <- as.integer(Dtau / Dt)
      Ntau_est <- as.integer(T_val / Dtau)
      NtNtau <- as.integer(T_val / Dt)

      incProgress(0.1, detail = "Parsing functions")
      f <- eval(parse(text = input$f_expr))
      h <- eval(parse(text = input$h_expr))
      df <- generate_derivative(f)

      incProgress(0.2, detail = "Constructing grid")
      x <- simData$x
      # Using global min/max; for different ranges consider per-dimension grids
      s_seq <- seq(min(x), max(x) + Ds, by = Ds)
      Ns <- length(s_seq)
      s_grid <- ExpandGrid(Dim, s_seq)

      incProgress(0.3, detail = "Computing matrices")
      D      <- generateD(Dim, Ns, Ds)
      B      <- computeB(s_grid, D, Dt, Ds, f, df, h)
      Lambda <- computeLambda(Dim, Ns, Dt, Ds)

      incProgress(0.3, detail = "Computing Iu")
      Iu <- wrap_outiu_function(s_grid, NtNtau, Ntau_est, Nt, Dim, simData$y, h, Lambda, B, Ns, NormalizedExp, DST_Solver)

      simData$Iu <- Iu
      incProgress(0.1, detail = "Completed")
    })

    cat("Iu estimation completed.\n")
  })

  # Dynamically generate the state trajectory tabs (one per dimension)
  output$state_plots <- renderUI({
    req(simData$x)
    Dim <- input$Dim
    tabs <- lapply(1:Dim, function(i) {
      tabPanel(title = paste0("Dimension ", i), value = i, plotOutput(paste0("state_plot_", i)))
    })
    do.call(tabsetPanel, c(list(id = "statePlots"), tabs))
  })

  # Dynamically generate the observation tabs (one per dimension)
  output$obs_plots <- renderUI({
    req(simData$y)
    Dim <- input$Dim
    tabs <- lapply(1:Dim, function(i) {
      tabPanel(title = paste0("Dimension ", i), value = i, plotOutput(paste0("obs_plot_", i)))
    })
    do.call(tabsetPanel, c(list(id = "obsPlots"), tabs))
  })

  # Render state trajectory plots
  observe({
    req(simData$x)
    Dim <- input$Dim
    for(i in 1:Dim) {
      local({
        idx <- i
        output[[paste0("state_plot_", idx)]] <- renderPlot({
          df_state <- data.frame(time = 1:nrow(simData$x), value = simData$x[, idx], type = "State")
          p <- ggplot() +
            geom_line(data = df_state, aes(x = time, y = value, color = type)) +
            labs(title = paste0("State Trajectory - Dimension ", idx),
                 x = "Time", y = paste0("Dimension ", idx)) +
            theme(
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA),
              panel.grid.major = element_line(colour = "lightgray"),
              panel.grid.minor = element_line(colour = "lightgray"),
              plot.background = element_blank(),
              legend.position= c(0.9, 0.9),  # 使用相对坐标将图例放置在右上角框线内
              legend.direction = "vertical",
              legend.background = element_rect(fill = "white", colour = "gray"),
              legend.box.background = element_rect(fill = "white", colour = "gray"),
              plot.title = element_text(hjust = 0.5, margin = margin(b = 10)),  # 调整标题边距
              legend.key.size = unit(0.5, "cm")  # 调整图例键大小
            ) +
            scale_y_continuous(labels = function(x) format(round(x, 1), nsmall = 1)) +
            scale_color_manual(name = NULL, values = c("State" = "#0099FF"))

          if (!is.null(simData$Iu)) {
            df_iu <- data.frame(time = 1:nrow(simData$Iu), value = simData$Iu[, idx], type = "YYAL Iu")
            p <- p +
              geom_line(data = df_iu, aes(x = time, y = value, color = type)) +
              scale_color_manual(name = NULL, values = c("State" = "#0099FF", "YYAL Iu" = "#FF3333"))
          }
          p
        })
      })
    }
  })

  # Render observation plots
  observe({
    req(simData$y)
    Dim <- input$Dim
    for(i in 1:Dim) {
      local({
        idx <- i
        output[[paste0("obs_plot_", idx)]] <- renderPlot({
          df_obs <- data.frame(time = 1:nrow(simData$y), value = simData$y[, idx], type = "Observation")
          ggplot(df_obs, aes(x = time, y = value, color = type)) +
            geom_point() +
            labs(title = paste0("Observation - Dimension ", idx),
                 x = "Time", y = paste0("Dimension ", idx)) +
            theme(
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA),
              panel.grid.major = element_line(colour = "lightgray"),
              panel.grid.minor = element_line(colour = "lightgray"),
              plot.background = element_blank(),
              legend.position= c(0.9, 0.9),  # 使用相对坐标将图例放置在右上角框线内
              legend.direction = "vertical",
              legend.background = element_rect(fill = "white", colour = "gray"),
              legend.box.background = element_rect(fill = "white", colour = "gray"),
              plot.title = element_text(hjust = 0.5, margin = margin(b = 10)),  # 调整标题边距
              legend.key.size = unit(0.5, "cm")  # 调整图例键大小
            ) +
            scale_y_continuous(labels = function(x) format(round(x, 1), nsmall = 1)) +
            scale_color_manual(name = NULL, values = c("Observation" = "#990000"))
        })
      })
    }
  })

  # Download handler to save the current plot as a PDF.
  # It checks which main tab is active ("State Trajectories" or "Observations")
  # and then which inner tab is selected.
  output$downloadPlot <- downloadHandler(
    filename = function() {
      if (input$mainTab == "State Trajectories") {
        paste0("state_plot_dimension_", input$statePlots, ".pdf")
      } else {
        paste0("observation_plot_dimension_", input$obsPlots, ".pdf")
      }
    },
    content = function(file) {
      if (input$mainTab == "State Trajectories") {
        idx <- as.numeric(input$statePlots)
        df_state <- data.frame(time = 1:nrow(simData$x), value = simData$x[, idx], type = "State")
        p <- ggplot() +
          geom_line(data = df_state, aes(x = time, y = value, color = type)) +
          labs(title = paste0("State Trajectory - Dimension ", idx),
               x = "Time", y = paste0("Dimension ", idx))+
          theme(
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA),
            panel.grid.major = element_line(colour = "lightgray"),
            panel.grid.minor = element_line(colour = "lightgray"),
            plot.background = element_blank(),
            legend.position = c(0.9, 0.9),  # 使用相对坐标将图例放置在右上角框线内
            legend.direction = "vertical",
            legend.background = element_rect(fill = "white", colour = "gray"),
            legend.box.background = element_rect(fill = "white", colour = "gray"),
            plot.title = element_text(hjust = 0.5, margin = margin(b = 10)),  # 调整标题边距
            legend.key.size = unit(0.5, "cm")  # 调整图例键大小
          ) +
          scale_y_continuous(labels = function(x) format(round(x, 1), nsmall = 1)) +
          scale_color_manual(name = NULL, values = c("State" = "#0099FF"))
        if (!is.null(simData$Iu)) {
          df_iu <- data.frame(time = 1:nrow(simData$Iu), value = simData$Iu[, idx], type ="YYAL Iu")
          p <- p +
            geom_line(data = df_iu, aes(x = time, y = value, color = type)) +
            scale_color_manual(name = NULL, values = c("State" = "#0099FF", "YYAL Iu" = "#FF3333"))
        }
        ggsave(file, plot = p, device = "pdf", width = 8, height = 6)
      } else if (input$mainTab == "Observations") {
        idx <- as.numeric(input$obsPlots)
        df_obs <- data.frame(time = 1:nrow(simData$y), value = simData$y[, idx], type = "Observation")
        p <- ggplot(df_obs, aes(x = time, y = value, color = type)) +
          geom_point() +
          labs(title = paste0("Observation - Dimension ", idx),
               x = "Time", y = paste0("Dimension ", idx)) +
          theme(
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA),
            panel.grid.major = element_line(colour = "lightgray"),
            panel.grid.minor = element_line(colour = "lightgray"),
            plot.background = element_blank(),
            legend.position = c(0.9, 0.9),  # 使用相对坐标将图例放置在右上角框线内
            legend.direction = "vertical",
            legend.background = element_rect(fill = "white", colour = "gray"),
            legend.box.background = element_rect(fill = "white", colour = "gray"),
            plot.title = element_text(hjust = 0.5, margin = margin(b = 10)),  # 调整标题边距
            legend.key.size = unit(0.5, "cm")  # 调整图例键大小
          ) +
          scale_y_continuous(labels = function(x) format(round(x, 1), nsmall = 1)) +
          scale_color_manual(name = NULL, values = c("Observation" = "#990000"))
        ggsave(file, plot = p, device = "pdf", width = 8, height = 6)
      }
    }
  )
}

shinyApp(ui, server)
