library(shiny)
library(bslib)
library(shinyjs)
library(DT)
library(plotly)
library(VEMIRT)
require(callr)
require(tibble)
require(openxlsx)
library(viridis)
library(shinyWidgets)
options(DT.options = list(scrollX = T))
source("exampledata/D2PL_rnd_gvem/shiny app D2PL_rnd_gvem all functions.R")
addResourcePath("static","www")
# input example: lambda 2-10 by 0.1 iteration = 50 MCMC iteration = 50 
# UI ----
ui <- page_navbar(
  window_title = 'VEMIRT',
  title = tags$a(img(src = 'static/logo.png', height="40px"), 'VEMIRT'),
  theme = bs_theme(bootswatch = 'litera', primary = "#007bff"),
  useShinyjs(),  # Enable shinyjs for sidebar toggling
  tags$head(
    tags$style(HTML("
    button.navbar-toggler,
    .navbar-collapse {
      display: none !important;
      visibility: hidden !important;
      width: 0 !important;
      height: 0 !important;
      overflow: hidden !important;
    }
  "))
  ),
  # Main Layout with Sidebar and Main Panel
  layout_sidebar(
    sidebar = sidebar(
      width = 300,  # Set sidebar width
      accordion(
        id = 'accord_sidebar',
        open = c('Method Selection', 'Data'),
        accordion_panel(
          'Method Selection',
          selectInput("method_selection", NULL,
                      choices = c("D2PL_pair_em" = "method1",
                                  "D2PL_rnd_gvem" = "method2"))
        ),
        
        accordion_panel(
          'Data',
          fileInput("fileUpload", "Upload a Data File (csv):", accept = ".csv"),
          fileInput("GroupUpload", "Upload a Group File (csv):", accept = ".csv")
        ),
        accordion_panel(
          'Parameters',
          
          # ========== For Method 1 ==========
          conditionalPanel(
            condition = "input.method_selection == 'method1'",
            
            tags$label(HTML("&lambda; Min:")),
            numericInput("lambda_min", NULL, value = 1, min = 0, step = 0.1),
            tags$label(HTML("&lambda; Max:")),
            numericInput("lambda_max", NULL, value = 1.5, min = 0, step = 0.1),
            tags$label(HTML("&lambda; By:")),
            numericInput("lambda_by", NULL, value = 0.1, min = 0.01, step = 0.01),
            
            tags$label(HTML("&tau; Min:")),
            numericInput("tau_min", NULL, value = 0.002, min = 0, step = 0.001),
            tags$label(HTML("&tau; Max:")),
            numericInput("tau_max", NULL, value = 0.01, min = 0, step = 0.001),
            tags$label(HTML("&tau; By:")),
            numericInput("tau_by", NULL, value = 0.002, min = 0.001, step = 0.001)
          ),
          
          # ========== For Method 2 ==========
          conditionalPanel(
            condition = "input.method_selection == 'method2'",
            
            tags$label(HTML("&lambda; Min:")),
            numericInput("lambda_min", NULL, value = 0, min = 0, step = 0.1),
            tags$label(HTML("&lambda; Max:")),
            numericInput("lambda_max", NULL, value = 20, min = 0, step = 0.1),
            tags$label(HTML("&lambda; By:")),
            numericInput("lambda_by", NULL, value = 0.1, min = 0.01, step = 0.01),
            
            tags$label("c-value:"),
            numericInput("cvalue", NULL, value = 0.04, min = 0, step = 0.001),
            
            tags$label("DIF Type:"),
            selectInput("dif_type", NULL, choices = c("UDIF" = "udif", "NON-UDIF" = "nudif")),
            
            tags$label("Impact Type:"),
            selectInput("impact_type", NULL, choices = c("Random" = "random", "Fixed" = "fixed")),
            
            tags$label("Iteration Criteria:"),
            numericInput("iter_criteria", NULL, value = 500, min = 1, step = 10),
            
            tags$label("Tau Criteria:"),
            numericInput("tau_criteria", NULL, value = 0.001, min = 0.0001, step = 0.0001),
            
            tags$label("rho_N:"),
            numericInput("rho_N", NULL, value = 0.001, min = 0, step = 0.0001),
            
            tags$label("rho_N2:"),
            numericInput("rho_N2", NULL, value = 0.1, min = 0, step = 0.01),
            
            tags$label("rho_Na:"),
            numericInput("rho_Na", NULL, value = 0.001, min = 0, step = 0.0001),
            
            tags$label("rho_Na2:"),
            numericInput("rho_Na2", NULL, value = 0.1, min = 0, step = 0.01),
            
            tags$head(
              tags$style(HTML("
      /* Make the selected text (on the button) black */
      #anchorDropdown .dropdown-toggle,
      #anchorDropdown .filter-option-inner-inner {
        color: black !important;
      }

      /* Make the dropdown menu items black */
      #anchorDropdown .dropdown-menu li a {
        color: black !important;
      }

      /* Optional: light blue background for the dropdown box */
      #anchorDropdown .dropdown-toggle {
        background-color: #e6f3ff !important;
        border: 1px solid #66b3ff;
        border-radius: 6px;
      }
    "))
            ),
            tags$label("Anchor Item:"),
            uiOutput("anchorDropdown"),
            
            
            # ===== Conditional: only show if dif_type == "udif" =====
            conditionalPanel(
              condition = "input.dif_type == 'udif'",
              
              checkboxInput("MCMC", label = "Use MCMC?", value = TRUE),
              tags$label("MCMC Iterations:"),
              numericInput("MCMC_iter", NULL, value = 50, min = 0, step = 10000),

            )
          )
        ),
        
        div(
          actionButton(
            "generate", "Run",
            style = "background-color:#99CCFF; border-color:#99CCFF; color:#000;"
          ),
            style = " text-align: center; margin-top: 15px;"
        ),
        div(style = "display:none;", textOutput("keepAlive"))
      )
    ),
    
    # Main Panel with Accordion Layout
    accordion(
      id = 'accord_output',
      
      # Uploaded Data Section
      accordion_panel(
        'Uploaded Data',
        navset_tab(
          nav_panel("Uploaded Data",
                    div(
                      DTOutput("dataTable"),
                      style = "margin-bottom:7px;"
                    )
          ),
          nav_panel("Uploaded Group",
                    div(DTOutput("groupTable"), style = "margin-bottom:7px;")
          )
        )
      ),
      
      # Visualization Section
      # Method 1: Show all visualizations
      conditionalPanel(
        condition = "input.method_selection == 'method1'",
        navset_tab(
          id = "viz_tabs",
          nav_panel("Group Heatmap", plotlyOutput("heatmapPlot")),
          nav_panel("SubGroup Mean", plotlyOutput("meanPlot")),
          nav_panel("Group Item Heatmap", plotlyOutput("GroupItemPlot")),
          nav_panel("Item ICC", uiOutput("iccPlots"))
        )
      ),
      
      # Method 2: Only show Heatmap
      # Method 2 outputs
      conditionalPanel(
        condition = "input.method_selection == 'method2'",
        
        # UDIF 
        conditionalPanel(
          condition = "input.dif_type == 'udif' && input.MCMC == true",
          navset_tab(
            id = "viz_tabs_method2_udif1",
            nav_panel("Group Heatmap", plotlyOutput("heatmapPlot2")),
            nav_panel("GIC Trend", plotlyOutput("gicPlots2")),
            nav_panel("Item ICC", uiOutput("iccPlots2")),
            nav_panel("Main Effect b", dataTableOutput("maineffectb0"))
          )
        ),
        conditionalPanel(
          condition = "input.dif_type == 'udif' && input.MCMC == false",
          navset_tab(
            id = "viz_tabs_method2_udif2",
            nav_panel("Group Heatmap", plotlyOutput("heatmapPlot3")),
            nav_panel("GIC Trend", plotlyOutput("gicPlots3")),
            nav_panel("Main Effect b", dataTableOutput("maineffectb1"))
          )
        ),
        
        # NUDIF 情况
        conditionalPanel(
          condition = "input.dif_type == 'nudif'",
          navset_tab(
            id = "viz_tabs_method2_nudif",
            nav_panel("Group Heatmap", plotlyOutput("heatmapPlot4")),
            nav_panel("GIC Trend", plotlyOutput("gicPlots4")),
            nav_panel("Main Effect a", dataTableOutput("maineffecta")),
            nav_panel("Main Effect b", dataTableOutput("maineffectb"))
          )
        )
      )
      
      
    )
  )
)

# server ----
server <- function(input, output, session) {
  
  options(shiny.sanitize.warnings = TRUE)
  options(shiny.autoreload = TRUE)
  safeRun <- function(expr) {
    tryCatch(
      expr,
      error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        message("Error: ", e$message)   # console
        shinyjs::enable("generate")
        shinyjs::html("generate", "Run")
        return(list())  
      }
    )
  }
  
  # Reactive value to store uploaded data
  uploaded_data <- reactiveVal(NULL)
  
  uploaded_group <- reactiveVal(NULL)
  # Reactive value for heatmap data
  heatmap_data <- reactiveVal(NULL)
  
  # Reactive value for heatmap data 2
  heatmap_data_2 <- reactiveVal(NULL)
  
  # Reactive value for heatmap data 3
  heatmap_data_3 <- reactiveVal(NULL)
  
  # Reactive value for heatmap data 4
  heatmap_data_4 <- reactiveVal(NULL)
  
  # Generate data function (uses uploaded data if available)
  generate_data <- function(uploaded_data = NULL, uploaded_group = NULL, Lambda, tau) {
    if (!is.null(uploaded_data) && !is.null(uploaded_group)) {
      # Run the provided D2PL_data process on uploaded data
      data = as.matrix(uploaded_data)
      group = as.integer(unlist(uploaded_group))
      dat <- D2PL_pair_em(data = data, group = group,
                                          Lambda0 = Lambda, #seq(1, 1.5, by = 0.1),
                                          Tau = tau, #c(Inf, seq(0.002, 0.01, by = 0.002)),
                                          verbose = TRUE)
      #fit <- with(D2PL_data, D2PL_pair_em(data, group, Lambda0 = seq(1, 1.5, by = 0.1), Tau = c(Inf, seq(0.002, 0.01, by = 0.002)), verbose = FALSE))
      
      fit <- dat 
      d <- fit$fit$d.a != 0 | fit$fit$d.b != 0
      dif <- rowSums(d, dims = 2)
      data <- dif + t(dif)
      rownames(data) <- paste0(" ", 1:length(unique(group)))
      colnames(data) <- paste0(" ", 1:length(unique(group)))
      return(list(data, coef(dat)$a, coef(dat)$b,coef(dat)$Sigma, coef(dat)$Mu, d))
    } else {
      # Generate random data if no uploaded data
      showNotification("No Data")
    }
  }
  
  generate_data_2 <- function(uploaded_data = NULL, uploaded_group = NULL,lambda,
    cvalue,dif_type,impact_type,anchor,MCMC,MCMC_iter,iter_criteria,tau_criteria ,
    rho_N,rho_N2,rho_Na,rho_Na2
  )
  {
    if (!is.null(uploaded_data) && !is.null(uploaded_group)) {
      # Run the provided D2PL_data process on uploaded data
      group = as.matrix(uploaded_group)
      showNotification("Running")
      dat <- D2PL_rnd_gvem(resp=uploaded_data, group_matrix = group, lambda,
                            cvalue,dif_type,impact_type,anchor,MCMC,MCMC_iter,iter_criteria,tau_criteria ,
                            rho_N,rho_N2,rho_Na,rho_Na2
      )      #fit <- with(D2PL_data, D2PL_pair_em(data, group, Lambda0 = seq(1, 1.5, by = 0.1), Tau = c(Inf, seq(0.002, 0.01, by = 0.002)), verbose = FALSE))
      # user choose
      # dat <- D2PL_rnd_gvem (uploaded_data, group, all_lambda = seq(0,10,length.out = 20), 
      #                              cvalue = 0.04, # user can modify
      #                              dif_type = "udif", #user input
      #                              impact_type = "random",
      #                              MCMC = TRUE, #defaulf false if non udif setup already
      #                              MCMC_iter = 1e2, # user can modify
      #                              iter_criteria = 5e2, #user can modify
      #                              tau_criteria = 1e-3, #user can modify
      #                              rho_N = 1e-3, #user can modify
      #                              rho_N2 = 0.1, #user can modify
      #                              rho_Na = 1e-3, #user can modify
      #                              rho_Na2 = 0.1 #user can modify
      #                              )
      
      fit <- dat
      data <- sign(fit$gvem_optimal_lambda$sig2_b_j)
      a <- colMeans(fit$MCMC_est$a_js)[,1]
      b <- colMeans(fit$MCMC_est$b_js)
      gic <- lapply(fit$gvem_all_lambda, function(x) x[["GIC2"]])
      main_eff_b <- fit$gvem_optimal_lambda$b_j
      return(list(data, a, b, gic, lambda, main_eff_b))
    } else {
      # Generate random data if no uploaded data
      showNotification("No Data")
    }
  }
  
  generate_data_3 <- function(uploaded_data = NULL, uploaded_group = NULL,lambda,
                              cvalue,dif_type,impact_type,anchor, MCMC,MCMC_iter,iter_criteria,tau_criteria ,
                              rho_N,rho_N2,rho_Na,rho_Na2
  )
  {
    if (!is.null(uploaded_data) && !is.null(uploaded_group)) {
      # Run the provided D2PL_data process on uploaded data
      group = as.matrix(uploaded_group)
      showNotification("Running")
      dat <- D2PL_rnd_gvem(resp=uploaded_data, group_matrix = group, lambda,
                            cvalue,dif_type,impact_type,anchor, FALSE,0,iter_criteria,tau_criteria ,
                            rho_N,rho_N2,rho_Na,rho_Na2
      )      #fit <- with(D2PL_data, D2PL_pair_em(data, group, Lambda0 = seq(1, 1.5, by = 0.1), Tau = c(Inf, seq(0.002, 0.01, by = 0.002)), verbose = FALSE))
      # user choose
      # dat <- D2PL_rnd_gvem (uploaded_data, group, all_lambda = seq(0,10,length.out = 20), 
      #                              cvalue = 0.04, # user can modify
      #                              dif_type = "udif", #user input
      #                              impact_type = "random",
      #                              MCMC = TRUE, #defaulf false if non udif setup already
      #                              MCMC_iter = 1e2, # user can modify
      #                              iter_criteria = 5e2, #user can modify
      #                              tau_criteria = 1e-3, #user can modify
      #                              rho_N = 1e-3, #user can modify
      #                              rho_N2 = 0.1, #user can modify
      #                              rho_Na = 1e-3, #user can modify
      #                              rho_Na2 = 0.1 #user can modify
      #                              )
      
      fit <- dat
      data <- sign(fit$gvem_optimal_lambda$sig2_b_j)
      gic <- lapply(fit$gvem_all_lambda, function(x) x[["GIC2"]])
      main_eff_b <- fit$gvem_optimal_lambda$b_j
      return(list(data, gic, lambda, main_eff_b))
    } else {
      # Generate random data if no uploaded data
      showNotification("No Data")
    }
  }
  generate_data_4 <- function(uploaded_data = NULL, uploaded_group = NULL,lambda,
                              cvalue,dif_type,impact_type,anchor,MCMC,MCMC_iter,iter_criteria,tau_criteria ,
                              rho_N,rho_N2,rho_Na,rho_Na2
  )
  {
    if (!is.null(uploaded_data) && !is.null(uploaded_group)) {
      # Run the provided D2PL_data process on uploaded data
      group = as.matrix(uploaded_group)
      showNotification("Running")
      cat(anchor)
      # dat <- D2PL_rnd_gvem (resp=uploaded_data, group_matrix = group, lambda,
      #                       cvalue,dif_type,impact_type,anchor = anchor,FALSE,0,iter_criteria,tau_criteria ,
      #                       rho_N,rho_N2,rho_Na,rho_Na2
      # )      #fit <- with(D2PL_data, D2PL_pair_em(data, group, Lambda0 = seq(1, 1.1, by = 0.1), Tau = c(Inf, seq(0.002, 0.002, by = 0.002)), verbose = FALSE))
      # # user choose
      uploaded_data <- read.csv("/Users/yijun.cheng/Downloads/shiny He/D2PL_rnd_resp_nudif.csv")
      group <- as.matrix(read.csv("/Users/yijun.cheng/Downloads/shiny He/D2PL_rnd_groupmat.csv"))
      dat <- D2PL_rnd_gvem (uploaded_data, group, all_lambda = seq(1,10,length.out = 10),
                                    cvalue = 0.04, # user can modify
                                    dif_type = "nudif", #user input
                                    impact_type = "fixed",
                                    anchor = c(19,20),
                                    MCMC = FALSE, #defaulf false if non udif setup already
                                    MCMC_iter = 1e3, # user can modify
                                    iter_criteria = 5e2, #user can modify
                                    tau_criteria = 1e-3, #user can modify
                                    rho_N = 1e-3, #user can modify
                                    rho_N2 = 0.1, #user can modify
                                    rho_Na = 1e-3, #user can modify
                                    rho_Na2 = 0.1 #user can modify
      )
      
      fit <- dat 
      
      data_a <- sign(fit$gvem_optimal_lambda$sig2_bar_a_j)
      data_b <- sign(fit$gvem_optimal_lambda$sig2_b_j)
      gic <- lapply(fit$gvem_all_lambda, function(x) x[["GIC2"]])
      main_eff_a <- fit$gvem_optimal_lambda$bar_a_j
      main_eff_b <- fit$gvem_optimal_lambda$b_j
      return(list(data_a,data_b, gic,lambda,main_eff_a, main_eff_b))
    } else {
      # Generate random data if no uploaded data
      showNotification("No Data")
    }
  }
  
  generate_list_from_matrix <- function(matrix) {
    n <- nrow(matrix)
    values <- seq_len(n)
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (matrix[i, j]) {
          #  matrix[i, j]  TRUE， i and j differ
          if (values[i] == values[j]) {
            values[j] <- max(values) + 1  # j new
          }
        } else {
          # matrix[i, j]  FALSE， i j same
          values[values == values[j]] <- values[i]
        }
      }
    }
    
    return(values)
  }
  # Handle file upload
  observeEvent(input$fileUpload, {
    req(input$fileUpload)  # Ensure a file is uploaded
    #data <- read.csv(input$fileUpload$datapath, row.names = 1)  # Read CSV file with row names
    data <- read.csv(input$fileUpload$datapath)
    uploaded_data(data)  # Store uploaded data as a matrix
    #runjs("$('#accord_output [data-value=\"Uploaded Data\"]').addClass('active').addClass('show');")
     
  })
  
  observe({
    data <- uploaded_data()
    req(data)
    
    output$anchorDropdown <- renderUI({
      pickerInput(
        inputId = "anchorItem",
        label = NULL,
        choices = names(data),
        selected = NULL,
        multiple = TRUE,
        options = pickerOptions(
          actionsBox = TRUE,
          liveSearch = TRUE,
          noneSelectedText = "Select columns as Anchor Items"
        )
      )
    })
  })
  
  
  # # --- Apply Selected Columns ---
  # observeEvent(input$applyCols, {
  #   data <- uploaded_data()
  #   req(data, input$selectedCols)
  #   
  #   filtered_data <- data[, input$selectedCols, drop = FALSE]
  #   uploaded_data(filtered_data)
  #   showNotification(
  #     paste("Applied column selection:",
  #           length(input$selectedCols), "columns kept"),
  #     type = "message"
  #   )
  # })
  
  observe({
    invalidateLater(60000, session)
    paste("Last keep-alive:", Sys.time())
  })
  # Handle file upload for the group file
  observeEvent(input$GroupUpload, {
    req(input$GroupUpload)  # Ensure a file is uploaded
    #group_data <- read.csv(input$GroupUpload$datapath, stringsAsFactors = FALSE)  # Read CSV file
    group_data <- read.csv(input$GroupUpload$datapath, stringsAsFactors = FALSE)
    uploaded_group(group_data)  # Update the reactive value
    #runjs("$('#accord_output [data-value=\"Uploaded Group\"]').addClass('show');")
    
    
  })
  # Handle default data generation or using uploaded data
  
  observeEvent(input$MCMC, {
    if (input$MCMC) {
      enable("MCMC_iter")
    } else {
      disable("MCMC_iter")
    }
  })
  observeEvent(input$generate, {
    gc()
    heatmap_data(NULL)
    heatmap_data_2(NULL)
    heatmap_data_3(NULL)
    heatmap_data_4(NULL)
    safeRun({
    # upload data
    data_file  <- input$fileUpload
    group_file <- input$GroupUpload
    
    if (input$method_selection == "method1") {
      shinyjs::html("generate", "Running...")
      shinyjs::disable("generate")
      # --- Method 1  ---
      lambda_min <- input$lambda_min
      lambda_max <- input$lambda_max
      lambda_by  <- input$lambda_by
      
      tau_min <- input$tau_min
      tau_max <- input$tau_max
      tau_by  <- input$tau_by
      
      Lambda0 <- seq(lambda_min, lambda_max, by = lambda_by)
      Tau     <- c(Inf, seq(tau_min, tau_max, by = tau_by))
      
      # generate data
      data <- generate_data(uploaded_data(), uploaded_group(), Lambda0, Tau)
      heatmap_data(data)  # update heatmap 
      #runjs("$('#accord_output [data-value=\"Uploaded Data\"]').removeClass('active');")
      #runjs("$('#accord_output [data-value=\"Uploaded Group\"]').removeClass('active');")
      #runjs("$('#accord_output [data-value=\"Visualizations\"]').addClass('active').addClass('show');")
      shinyjs::enable("generate")
      shinyjs::html("generate", "Run")
      
    } else if (input$method_selection == "method2") {
      # --- Method 2  ---
      
      lambda_min <- input$lambda_min
      lambda_max <- input$lambda_max
      lambda_by  <- input$lambda_by
      
      all_lambda     <- seq(lambda_min, lambda_max, by = lambda_by)
      cvalue         <- input$cvalue
      dif_type       <- input$dif_type
      impact_type    <- input$impact_type
      MCMC           <- ifelse(dif_type == "udif", input$MCMC, FALSE)  # when udif: MCMC is TRUE
      MCMC_iter      <- input$MCMC_iter
      iter_criteria  <- input$iter_criteria
      tau_criteria   <- input$tau_criteria
      rho_N          <- input$rho_N
      rho_N2         <- input$rho_N2
      rho_Na         <- input$rho_Na
      rho_Na2        <- input$rho_Na2
      anchors <- input$anchorItem
      data <- uploaded_data()
      anchor_indices <- which(names(data) %in% anchors)
      print(class(anchor_indices))
      # conduct method2 
      if(input$dif_type == "udif"){
      shinyjs::html("generate", "Running...")
      shinyjs::disable("generate")
      if (input$MCMC) {
        data <- generate_data_2(
          uploaded_data(),
          uploaded_group(),
          all_lambda,cvalue,
          dif_type,impact_type,
          anchor = anchor_indices,
          MCMC  ,
          MCMC_iter ,
          iter_criteria ,
          tau_criteria ,
          rho_N,
          rho_N2 ,
          rho_Na  ,
          rho_Na2
        )
        heatmap_data_2(data)
        } else {
          data <- generate_data_3(
            uploaded_data(),
            uploaded_group(),
            all_lambda,cvalue,
            dif_type,impact_type,
            anchor = anchor_indices,
            FALSE,
            0,
            iter_criteria ,
            tau_criteria ,
            rho_N,
            rho_N2 ,
            rho_Na  ,
            rho_Na2
          )
          heatmap_data_3(data)
      }
   

      #runjs("$('#accord_output [data-value=\"Uploaded Data\"]').removeClass('active');")
      #runjs("$('#accord_output [data-value=\"Uploaded Group\"]').removeClass('active');")
      #runjs("$('#accord_output [data-value=\"Visualizations\"]').addClass('active').addClass('show');")
      shinyjs::enable("generate")
      shinyjs::html("generate", "Run")
      
      }else{
        shinyjs::html("generate", "Running...")
        shinyjs::disable("generate")
        data <- generate_data_4(
          uploaded_data(),
          uploaded_group(),
          all_lambda,cvalue,
          dif_type,impact_type,
          anchor = anchor_indices,
          FALSE  ,
          0 ,
          iter_criteria ,
          tau_criteria ,
          rho_N,
          rho_N2 ,
          rho_Na  ,
          rho_Na2
        )
        heatmap_data_4(data)  
        #runjs("$('#accord_output [data-value=\"Uploaded Data\"]').removeClass('active');")
        #runjs("$('#accord_output [data-value=\"Uploaded Group\"]').removeClass('active');")
        #runjs("$('#accord_output [data-value=\"Visualizations\"]').addClass('active').addClass('show');")
        shinyjs::enable("generate")
        shinyjs::html("generate", "Run")
        
      }
        # 
    }
    
     })
  })
  
  # Render data table (reflects uploaded data only)
  output$dataTable <- renderDT({
    data <- uploaded_data()  # Get the uploaded data
    if (is.null(data)) {
      datatable(data.frame(), options = list(pageLength = 10),
                caption = "No data available. Upload a file to view the data.")
    } else {
      datatable(data,  options = list(
        scrollX = "200px",  
        scrollY = "350px", 
        paging = FALSE,  
        searching = TRUE,  
        fixedColumns = list(leftColumns = 1)  
      ))
    }
  })
  # Render the uploaded group file in a data table
  output$groupTable <- renderDT({
    data <- uploaded_group()  # Get the uploaded group file data
    if (is.null(data)) {
      datatable(data.frame(), options = list(pageLength = 10), 
                caption = "No data available. Upload a group file to view the data.")
    } else {
      datatable(data, options = list(
        scrollX = "200px",  
        scrollY = "350px",  
        paging = FALSE,  
        searching = TRUE,  
        fixedColumns = list(leftColumns = 1)
      ))
    }
  })
  # Render heatmap
  output$heatmapPlot <- renderPlotly({
    data <- heatmap_data() # Get the heatmap data
    req(data)  # Ensure heatmap data is not NULL
    req(!is.null(data) && length(data) > 0)
    data <- data[[1]]
    heatmap_df <- expand.grid(
      Row = rownames(data),
      Col = colnames(data)
    )
    heatmap_df$Dif <- as.vector(data)  # Flatten the matrix into a single column for plotting
    
    # Generate heatmap with hover information
    plot_ly(
      data = heatmap_df,
      x = ~Row,
      y = ~Col,
      z = ~Dif,
      type = "heatmap",
      colors = c("blue", "yellow"),
      hoverinfo = "text",
      text = ~paste(
        "Group: ", Row, "<br>",
        "Group: ", Col, "<br>",
        "DIF: ", round(Dif, 2)
      )
    ) %>%
      layout(
        title = "Interactive Heatmap DIF",
        xaxis = list(title = "Group"),
        yaxis = list(title = "Group"),
        hoverlabel = list(bgcolor = "white", font = list(color = "black"))
      ) %>% config(
        modeBarButtonsToRemove = list(
          "pan2d", "select2d", "lasso2d", "hoverClosestCartesian", 
          "hoverCompareCartesian", "autoScale2d", 
          "resetScale2d", "toggleSpikelines",  # Remove reset zoom
          "zoom2d"  # Remove zoom tool
        ),
        displaylogo = FALSE  # Remove Plotly logo
      )
  })
  
  output$heatmapPlot2 <- renderPlotly({
    data <- heatmap_data_2() # Get the heatmap data
    req(data)  # Ensure heatmap data is not NULL
    validate(need(!is.null(data), "No data available."))
    req(!is.null(data) && length(data) > 0)
    data <- data[[1]]
    colors <- viridis(2, alpha = 0.7)
    Binary = ifelse(data, "TRUE", "FALSE")
    position = 1:length(data)
    heatmap_df <- data.frame(
      position = position,
      Binary = Binary,
      HoverInfo = paste("Item:", position, "<br>DIF:", Binary),
      Color = colors[data+1]
    )
    p <- ggplot(heatmap_df, aes(x = 1, y = (position))) +
      geom_tile(aes(fill = Color, text = HoverInfo), color = "white") +  # Use Color column for fill
      scale_fill_identity() +  # Directly use colors from Color column
      labs(
        title = "Interactive Heatmap DIF",
        x = "DIF",
        y = "Item"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.text.x=element_blank(),
        legend.position = "none"  # Remove legend
      )
    # Convert ggplot to plotly
    ggplotly(p, tooltip = c("text")) 
  })
  
  output$heatmapPlot3 <- renderPlotly({
    data <- heatmap_data_3() # Get the heatmap data
    req(data)  # Ensure heatmap data is not NULL
    validate(need(!is.null(data), "No data available."))
    req(!is.null(data) && length(data) > 0)
    data <- data[[1]]
    colors <- viridis(2, alpha = 0.7)
    Binary = ifelse(data, "TRUE", "FALSE")
    position = 1:length(data)
    heatmap_df <- data.frame(
      position = position,
      Binary = Binary,
      HoverInfo = paste("Item:", position, "<br>DIF:", Binary),
      Color = colors[data+1]
    )
    p <- ggplot(heatmap_df, aes(x = 1, y = (position))) +
      geom_tile(aes(fill = Color, text = HoverInfo), color = "white") +  # Use Color column for fill
      scale_fill_identity() +  # Directly use colors from Color column
      labs(
        title = "Interactive Heatmap DIF",
        x = "DIF",
        y = "Item"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.text.x=element_blank(),
        legend.position = "none"  # Remove legend
      )
    # Convert ggplot to plotly
    ggplotly(p, tooltip = c("text")) 
  })
  
  
  output$heatmapPlot4 <- renderPlotly({
    data <- heatmap_data_4() # Get the heatmap data
    req(data)  # Ensure heatmap data is not NULL
    req(!is.null(data) && length(data) > 0)
    data1 <- data[[1]]
    data2 <- data[[2]]
    colors <- viridis(2, alpha = 0.7)
    Binary1 = ifelse(data1, "TRUE", "FALSE")
    Binary2 = ifelse(data2, "TRUE", "FALSE")
    position1 = 1:length(data1)
    position2 = 1:length(data2)
    heatmap_df1 <- data.frame(
      position = position1,
      Binary = Binary1,
      HoverInfo = paste("Item:", position1, "<br>DIF on a:", Binary1),
      Color = colors[data1+1]
    )
    heatmap_df2 <- data.frame(
      position = position2,
      Binary = Binary2,
      HoverInfo = paste("Item:", position2, "<br>DIF on b:", Binary2),
      Color = colors[data2+1]
    )
    
    heatmap_df1$Source <- "Interactive Heatmap DIF on a"
    heatmap_df2$Source <- "Interactive Heatmap DIF on b"
    df_all <- rbind(heatmap_df1, heatmap_df2)
    
    p <- ggplot(df_all, aes(x = 1, y = position)) +
      geom_tile(aes(fill = Color, text = HoverInfo), color = "white") +
      scale_fill_identity() +
      facet_wrap(~Source) +
      labs(title = "Interactive Heatmap DIF", x = "", y = "Item") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        legend.position = "none"
      )
    # Convert ggplot to plotly
    ggplotly(p, tooltip = c("text")) 
  })
  
  # Generate ICC
  output$iccPlots <- renderUI({
    data <- heatmap_data()  # Get the heatmap data
    req(data)  # Ensure heatmap data is not NULL
    req(!is.null(data) && length(data) > 0)
    
    # Simulated data for a_values and b_values (replace with your real data)
    a_values <- data[[2]]  # Discrimination parameters (rows: groups, cols: items)
    b_values <- data[[3]]  # Difficulty parameters (rows: groups, cols: items)
    
    n_items <- ncol(a_values)
    theta <- seq(-3, 3, length.out = 100)  # Define ability levels
    
    # Use a predefined color palette with transparency
    transparent_colors <- viridis(nrow(a_values), alpha = 0.3)
    # Generate a Plotly output for each item
    plot_outputs <- lapply(seq_len(n_items), function(i) {
      plotname <- paste0("plot_", i)
      output[[plotname]] <- renderPlotly({
        plot_data <- data.frame()
        
        p <- plot_ly()  # Initialize an empty plot
        
        for (j in seq_len(nrow(a_values))) {  # Loop over groups
          a <- a_values[j, i]
          b <- b_values[j, i]
          prob <- 1 / (1 + exp(-a * (theta - b)))
          
          curve_data <- data.frame(
            Theta = as.vector(theta),  # Ensure no row names are carried over
            Probability = as.vector(prob),
            Curve = paste0("Group ", j),
            HoverInfo = paste(
              "Group:", j,
              "<br>a=", round(a, 2),
              "<br>b=", round(b, 2)
            )
          )
          
          # Add trace for each group
          p <- p %>%
            add_trace(
              data = curve_data,
              x = ~Theta,
              y = ~Probability,
              type = 'scatter',
              mode = 'lines',
              line = list(color = transparent_colors[j], width = 2.5),  # Transparent color
              name = paste0("Group ", j),  # Legend entry
              text = ~HoverInfo,
              hoverinfo = "text"
            )
        }
        
        p %>%
          layout(
            title = paste("Item", i, "Characteristic Curves (ICCs)"),
            xaxis = list(title = expression(theta)),
            yaxis = list(title = "Probability of Correct Response"),
            legend = list(title = list(text = "Groups"))
          ) %>%
          config(
            modeBarButtonsToRemove = list(
              "pan2d", "select2d", "lasso2d", "hoverClosestCartesian", 
              "hoverCompareCartesian", "autoScale2d", 
              "resetScale2d", "toggleSpikelines",  # Remove reset zoom
              "zoom2d"  # Remove zoom tool
            ),
            displaylogo = FALSE  # Remove Plotly logo
          )
      })
      
      div(
        plotlyOutput(plotname, height = "300px"),
        style = "margin-bottom: 30px;"
      )
    })
    
    do.call(tagList, plot_outputs)
  })
  
  output$iccPlots2 <- renderUI({
    data <- heatmap_data_2()  # Get the heatmap data
    req(data)  # Ensure heatmap data is not NULL
    req(!is.null(data) && length(data) > 0)
    
    # Simulated data for a_values and b_values (replace with your real data)
    a_values <- data[[2]]  # Discrimination parameters (items)
    b_values <- t(data[[3]])  # Difficulty parameters (rows: groups, cols: items)
    n_items <- length(a_values)
    theta <- seq(-3, 3, length.out = 100)  # Define ability levels
    
    # Use a predefined color palette with transparency
    transparent_colors <- viridis(nrow(b_values), alpha = 0.3)
    # Generate a Plotly output for each item
    plot_outputs_2 <- lapply(seq_len(n_items), function(i) {
      plotname_2 <- paste0("plot_", i, "_section2") 
      output[[plotname_2]] <- renderPlotly({
        plot_data <- data.frame()
        
        p <- plot_ly()  # Initialize an empty plot
        
        for (j in seq_len(length(a_values))) {  # Loop over groups
          a <- a_values[i]
          b <- b_values[j, i]
          prob <- 1 / (1 + exp(-a * (theta - b)))
          
          curve_data <- data.frame(
            Theta = as.vector(theta),  # Ensure no row names are carried over
            Probability = as.vector(prob),
            Curve = paste0("Group ", j),
            HoverInfo = paste(
              "Group:", j,
              "<br>a=", round(a, 2),
              "<br>b=", round(b, 2)
            )
          )
          
          # Add trace for each group
          p <- p %>%
            add_trace(
              data = curve_data,
              x = ~Theta,
              y = ~Probability,
              type = 'scatter',
              mode = 'lines',
              line = list(color = transparent_colors[j], width = 2.5),  # Transparent color
              name = paste0("Group:", j, ";b:",round(b, 2)),  # Legend entry
              text = ~HoverInfo,
              hoverinfo = "text"
            )
        }
        
        p %>%
          layout(
            title = paste("Item", i, "Characteristic Curves (ICCs)"),
            xaxis = list(title = expression(theta)),
            yaxis = list(title = "Probability of Correct Response"),
            legend = list(title = list(text = "Groups"))
          ) %>%
          config(
            modeBarButtonsToRemove = list(
              "pan2d", "select2d", "lasso2d", "hoverClosestCartesian", 
              "hoverCompareCartesian", "autoScale2d", 
              "resetScale2d", "toggleSpikelines",  # Remove reset zoom
              "zoom2d"  # Remove zoom tool
            ),
            displaylogo = FALSE  # Remove Plotly logo
          )
      })
      
      div(
        plotlyOutput(plotname_2, height = "300px"),
        style = "margin-bottom: 30px;"
      )
    })
    
    do.call(tagList, plot_outputs_2)
  })
  
  output$GroupItemPlot <- renderPlotly({
    # Create ggplot heatmap
    data <- heatmap_data()  # Get the heatmap data
    req(data)  # Ensure heatmap data is not NULL
    req(!is.null(data) && length(data) > 0)
    
    # Simulated data for a_values and b_values (replace with your real data)
    a_values <- round(data[[2]],2)  # Discrimination parameters (rows: groups, cols: items)
    b_values <- round(data[[3]],2) # Difficulty parameters (rows: groups, cols: items)
    d <- data[[6]]
    # Prepare a data frame for plotting
    
    n_colors <- nrow(a_values)
    colors <- viridis(n_colors, alpha = 0.7)
    
    heatmap_data <- data.frame()
    for (i in seq_len(ncol(a_values))) {  # Loop through each item
      # Combine a and b for the current item
      index <- generate_list_from_matrix(d[,,i])
      unique_combinations <- unique(index)
      # Generate a color palette for unique combinations
      color_map <- setNames(colors, unique_combinations)
      # Add data for the current item
      heatmap_data <- rbind(
        heatmap_data,
        data.frame(
          Item = rep(i, times = nrow(a_values)),  # Rows represent items
          Group = seq_len(nrow(a_values)),  # Columns represent groups
          Value = index,
          HoverInfo = paste("Item:", i, 
                            "<br>Group:", seq_len(nrow(a_values)), 
                            "<br>a:", round(a_values[seq_len(nrow(a_values)),i],2),
                            "<br>b:", round(b_values[seq_len(nrow(a_values)),i],2)),
          Color = color_map[index]
        )
      )
    }
    heatmap_data <- heatmap_data %>%
      mutate(Item = factor(Item, levels = rev(sort(unique(Item)))))  # Sort Item in descending order
    
    p <- ggplot(heatmap_data, aes(x = factor(Group), y = factor(Item))) +
      geom_tile(aes(fill = Color, text = HoverInfo), color = "white") +  # Use Color column for fill
      scale_fill_identity() +  # Directly use colors from Color column
      labs(
        title = "Item vs Group",
        x = "Group",
        y = "Item"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text( hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        legend.position = "none"  # Remove legend
      )
    
    # Convert ggplot to plotly
    ggplotly(p, tooltip = c("text")) 
  })
  output$meanPlot <- renderPlotly({
    # Create ggplot heatmap
    data <- heatmap_data()  # Get the heatmap data
    req(data)  # Ensure heatmap data is not NULL
    req(!is.null(data) && length(data) > 0)
    
    # Simulated data for a_values and b_values (replace with your real data)
    sigma_values <- data[[4]]  # Discrimination parameters (rows: groups, cols: items)
    mu_values <- data[[5]]  # Difficulty parameters (rows: groups, cols: items)
    # Define x range for the plot
    numCurves <- length(sigma_values)
    x <- seq(-20, 20, length.out = 1000)  # Define x range
    
    # Initialize plotly object
    p <- plot_ly()
    color <- viridis(numCurves, alpha = 0.7)
    # Add each normal curve to the plot
    for (i in 1:numCurves) {
      if (!is.null(mu_values[i]) && !is.null(sigma_values[i])) {
        y <- dnorm(x, mean = mu_values[i], sd = sigma_values[i])  # Calculate normal density
        p <- p %>%
          add_trace(
            x = x,
            y = y,
            type = "scatter",
            mode = "lines",
            name = paste(i),
            line = list(color =  color[i], width = 2)
          )
      }
    }
    
    # Customize the layout
    p <- p %>%
      layout(
        title = "Dynamic Normal Curves",
        xaxis = list(title = "X"),
        yaxis = list(title = "Density"),
        legend = list(title = list(text = "Groups"))
      )%>% config(
        modeBarButtonsToRemove = list(
          "pan2d", "select2d", "lasso2d", "hoverClosestCartesian", 
          "hoverCompareCartesian", "autoScale2d", 
          "resetScale2d", "toggleSpikelines",  # Remove reset zoom
          "zoom2d"  # Remove zoom tool
        ),
        displaylogo = FALSE  # Remove Plotly logo
      )
    
    p
  })
  output$gicPlots2 <- renderPlotly({
    data <- heatmap_data_2()  # Get the heatmap data
    req(data)  # Ensure heatmap data is not NULL
    req(!is.null(data) && length(data) > 0)
    plot_ly(
      x = ~data[[5]], 
      y = ~data[[4]], 
      type = 'scatter', 
      mode = 'lines+markers',
      hoverinfo = "x+y",  # x y
      hoverlabel = list(
        bgcolor = "white",  # 
        bordercolor = "black",  #
        font = list(color = "black") 
      )
    ) %>%
      layout(
        title = "GIC Trend with Lambda",
        xaxis = list(title = "Lambda"),
        yaxis = list(title = "GIC"),
        hovermode = "closest",  #  hover 
        plot_bgcolor = "white",  # 
        paper_bgcolor = "white"  # 
      ) %>%
      config(
        modeBarButtonsToRemove = list(
          "pan2d", "select2d", "lasso2d", "hoverClosestCartesian",
          "hoverCompareCartesian", "autoScale2d",
          "resetScale2d", "toggleSpikelines",
          "zoom2d"
        ),
        displaylogo = FALSE
      )
  })
  output$gicPlots3 <- renderPlotly({
    data <- heatmap_data_3()  # Get the heatmap data
    req(data)  # Ensure heatmap data is not NULL
    req(!is.null(data) && length(data) > 0)
    
    df3 <- data.frame(lambda = data[[3]], gic = data[[2]])
    plot_ly(
      x = ~data[[3]], 
      y = ~data[[2]], 
      type = 'scatter', 
      mode = 'lines+markers',
      hoverinfo = "x+y",  # x y
      hoverlabel = list(
        bgcolor = "white",  # 
        bordercolor = "black",  #
        font = list(color = "black") 
      )
    ) %>%
      layout(
        title = "GIC Trend with Lambda",
        xaxis = list(title = "Lambda"),
        yaxis = list(title = "GIC"),
        hovermode = "closest",  #  hover 
        plot_bgcolor = "white",  # 
        paper_bgcolor = "white"  # 
      ) %>%
      config(
        modeBarButtonsToRemove = list(
          "pan2d", "select2d", "lasso2d", "hoverClosestCartesian",
          "hoverCompareCartesian", "autoScale2d",
          "resetScale2d", "toggleSpikelines",
          "zoom2d"
        ),
        displaylogo = FALSE
      )
  })
  output$gicPlots4 <- renderPlotly({
    data <- heatmap_data_4()  # Get the heatmap data
    req(data)  # Ensure heatmap data is not NULL
    req(!is.null(data) && length(data) > 0)
    
    df3 <- data.frame(lambda = data[[4]], gic = data[[3]])
    plot_ly(
      x = ~data[[4]], 
      y = ~data[[3]], 
      type = 'scatter', 
      mode = 'lines+markers',
      hoverinfo = "x+y",  # x y
      hoverlabel = list(
        bgcolor = "white",  # 
        bordercolor = "black",  #
        font = list(color = "black") 
      )
    ) %>%
      layout(
        title = "GIC Trend with Lambda",
        xaxis = list(title = "Lambda"),
        yaxis = list(title = "GIC"),
        hovermode = "closest",  #  hover 
        plot_bgcolor = "white",  # 
        paper_bgcolor = "white"  # 
      ) %>%
      config(
        modeBarButtonsToRemove = list(
          "pan2d", "select2d", "lasso2d", "hoverClosestCartesian",
          "hoverCompareCartesian", "autoScale2d",
          "resetScale2d", "toggleSpikelines",
          "zoom2d"
        ),
        displaylogo = FALSE
      )
  })
  output$maineffecta <- renderDT({
    data <- heatmap_data_4() # Get the heatmap data
    req(data)  # Ensure heatmap data is not NULL
    req(!is.null(data) && length(data) > 0)
    data_effect_a <- round(data[[5]],3)
    data_effect_a <- cbind(Row = seq_len(nrow(data_effect_a)), data_effect_a)
    new_names <- c("Row", paste0("G", seq_len(ncol(data_effect_a) - 1)))
    colnames(data_effect_a) <- new_names
    datatable(data_effect_a, options = list(pageLength = 20))
  })
  output$maineffectb <- renderDT({
    data <- heatmap_data_4() # Get the heatmap data
    req(data)  # Ensure heatmap data is not NULL
    req(!is.null(data) && length(data) > 0)
    data_effect_b <- round(data[[6]],3)
    data_effect_b <- cbind(Row = seq_len(nrow(data_effect_b)), data_effect_b)
    new_names <- c("Row", paste0("G", seq_len(ncol(data_effect_b) - 1)))
    colnames(data_effect_b) <- new_names
    datatable(data_effect_b, options = list(pageLength = 20))
  })
  output$maineffectb0 <- renderDT({
    data <- heatmap_data_2() # Get the heatmap data
    req(data)  # Ensure heatmap data is not NULL
    req(!is.null(data) && length(data) > 0)
    data_effect_b <- round(data[[6]],3)
    data_effect_b <- cbind(Row = seq_len(nrow(data_effect_b)), data_effect_b)
    new_names <- c("Row", paste0("G", seq_len(ncol(data_effect_b) - 1)))
    colnames(data_effect_b) <- new_names
    datatable(data_effect_b, options = list(pageLength = 20))
  })
  output$maineffectb1 <- renderDT({
    data <- heatmap_data_3() # Get the heatmap data
    req(data)  # Ensure heatmap data is not NULL
    req(!is.null(data) && length(data) > 0)
    data_effect_b <- round(data[[4]],3)
    data_effect_b <- cbind(Row = seq_len(nrow(data_effect_b)), data_effect_b)
    new_names <- c("Row", paste0("G", seq_len(ncol(data_effect_b) - 1)))
    colnames(data_effect_b) <- new_names
    datatable(data_effect_b, options = list(pageLength = 20))
  })
}

# Run Shiny App
shinyApp(ui, server)
