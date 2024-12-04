# Written by Weicong Lyu
require(shiny)
require(shinyjs)
require(bslib)
require(data.table)
require(DT)
require(callr)
require(tibble)
require(openxlsx)
require(VEMIRT)

options(DT.options = list(scrollX = T))

ui <- page_navbar(
  window_title = 'VEMIRT',
  title = tags$a(img(src = 'logo.png', height = '45px'), 'VEMIRT'),
  id = 'navbar',
  tags$style('
    .navbar { z-index: 99999 }
    td.dt-right {
      font-variant-numeric: tabular-nums;
    }
    input[type = "number"] {
      text-align: center;
      height: 30px;
      -moz-appearance: textfield;
    }
    input[type = "number"]::-webkit-inner-spin-button, input[type = "number"]::-webkit-outer-spin-button {
      -webkit-appearance: none;
      margin: 0;
    }
    input[type = "number"]::placeholder {
      font-style: italic;
    }
    .compact .shiny-input-container {
      margin-bottom: .5rem;
    }
  '),
  tags$script('
    $(document).on("shiny:connected", function() {
      $("#b_mu").attr("placeholder", "μ");
      $("#b_sigma2").attr("placeholder", "σ²");
      $("#c_alpha").attr("placeholder", "α");
      $("#c_beta").attr("placeholder", "β");
    });
  '),
  theme = bs_theme(bootswatch = 'litera'),
  useShinyjs(),
  nav_menu(
    'EFA',
    nav_panel('2PL', value = 1),
    nav_panel('3PL', value = 2)
  ),
  nav_menu(
    'CFA',
    nav_panel('2PL', value = 11),
    nav_panel('3PL', value = 12)
  ),
  nav_menu(
    'DIF',
    nav_panel('2PL', value = 21),
    #nav_panel('3PL', value = 22)
  ),
  nav_spacer(),
  nav_item(tags$a(icon('home'), href = 'https://MAP-LAB-UW.github.io/VEMIRT/', target = '_blank')),
  nav_item(tags$a(icon('github'), href = 'https://github.com/MAP-LAB-UW/VEMIRT', target = '_blank')),
  header = page_sidebar(
    sidebar = sidebar(
      width = 300,
      accordion(
        id = 'accord_sidebar',
        open = c('Data', 'Model'),
        accordion_panel(
          'Data',
          fileInput('in_data', 'Responses:'),
          conditionalPanel('input.navbar >= 10', fileInput('in_model', 'Loading indicators:')),
          conditionalPanel('input.navbar >= 20', fileInput('in_group', 'Group indicators:'))
        ),
        accordion_panel(
          'Model',
          selectInput('algorithm', 'Algorithm:', NULL),
          conditionalPanel(
            'input.navbar < 10',
            div(
              class = 'compact',
              style = 'display: inline-flex; align-items: stretch',
              tags$label(HTML('Number of factors:&nbsp;')),
              numericInput('nfactor', NULL, 3, min = 1, step = 1, width = '75px')
            ),
            conditionalPanel(
              'input.algorithm == "rotation"',
              selectInput('rotation', 'Post-hoc rotation method:', c('Promax', 'CF-Quartimax' = 'cfQ'), 'Promax')
            ),
            conditionalPanel(
              'input.algorithm != "rotation"',
              selectInput('constrain_type', 'Constraint:', c('Identity matrix' = 'C1', 'Lower triangular matrix' = 'C2'), 'C1'),
              uiOutput('constrain_items')
            )
          ),
          conditionalPanel(
            'input.navbar % 10 == 2',
            tags$label('Priors:'),
            div(
              class = 'compact',
              style = 'display: inline-flex; align-items: stretch',
              span(HTML('<i>b</i> ~ Normal(')),
              numericInput('b_mu', NULL, 0, step = 0.1, width = '75px'),
              span(', '),
              numericInput('b_sigma2', NULL, 4, min = 0, step = 1, width = '75px'),
              span(')')
            ),
            div(
              style = 'display: inline-flex; align-items: stretch',
              span(HTML('<i>c</i> ~ Beta(')),
              numericInput('c_alpha', NULL, 10, min = 0, step = 0.1, width = '75px'),
              span(', '),
              numericInput('c_beta', NULL, 40, min = 0, step = 0.1, width = '75px'),
              span(')')
            )
          ),
          conditionalPanel(
            'input.navbar == 21 && (input.algorithm == "EM" || input.algorithm == "EMM")',
            div(
              class = 'compact',
              style = 'display: inline-flex; align-items: stretch',
              tags$label(HTML('Accuracy level:&nbsp;')),
              numericInput('level', NULL, 10, 1, 25, 1, width = '75px')
            )
          ),
          conditionalPanel(
            'input.navbar >= 20',
            tags$label('Tuning parameters:'),
            div(
              style = 'display: flex; align-items: stretch',
              div(
                style = 'display: inline-flex; white-space: nowrap',
                tags$label(HTML('<i>λ</i><sub>0</sub> =&nbsp;'))
              ),
              div(
                style = 'flex: 0; text-align: right',
                div(
                  class = 'compact',
                  style = 'display: inline-flex',
                  numericInput('lambda_from', NULL, 0, min = 0, step = 0.1, width = '90px'),
                  span(HTML('&nbsp;to&nbsp;')),
                  numericInput('lambda_to', NULL, 1, min = 0, step = 0.1, width = '90px')
                ),
                div(
                  style = 'display: inline-flex',
                  span(HTML('by&nbsp;')),
                  numericInput('lambda_by', NULL, 0.2, min = 0, step = 0.01, width = '90px')
                )
              )
            )
          ),
          actionButton('in_run', 'Run')
        ),
        accordion_panel(
          'Postestimation',
          div(
            id = 'panel_ic',
            class = 'compact',
            style = 'display: none',
            selectizeInput('ic', 'Information Criterion:',
              list("AIC (<i>l</i><sub>0</sub>⋅2)" = 'AIC',
                   "BIC (<i>l</i><sub>0</sub>⋅log<i>N</i>)" = 'BIC',
                   "GIC (<i>l</i><sub>0</sub>⋅log<i>N</i>⋅<i>c</i>loglog<i>N</i>)" = 'GIC'),
              options = list(render = I('{
                item: function(item, escape) {
                  return "<div>" + item.label + "</div>"
                },
                option: function(item, escape) {
                  return "<div>" + item.label + "</div>"
                }
              }'))
            ),
            conditionalPanel(
              'input.ic == "GIC"',
              style = 'display: inline-flex; align-items: stretch',
              tags$label(HTML('<i>c</i> =&nbsp;')),
              numericInput('gic_c', NULL, 1, min = 0, step = 0.1, width = '90px')
            ),
            div(
              style = 'margin-bottom: 1rem',
              tags$label(HTML('<i>λ</i><sub>0</sub><sup>*</sup> =&nbsp;')),
              textOutput('out_lambda0', inline = T)
            )
          ),
          downloadButton('download')
        )
      )
    ),
    accordion(
      id = 'accord_output',
      accordion_panel(
        'Data',
        navset_tab(
          nav_panel('Responses', DTOutput('out_data')),
          nav_panel('Loading', DTOutput('out_model')),
          nav_panel('Group', DTOutput('out_group')),
        )
      ),
      accordion_panel(
        'Parameter Estimates',
        navset_tab(
          nav_panel('Slope', DTOutput('out_a')),
          nav_panel('DIF in Slope', DTOutput('out_gamma')),
          nav_panel('Intercept', DTOutput('out_b')),
          nav_panel('DIF in Intercept', DTOutput('out_beta')),
          nav_panel('Guessing', DTOutput('out_c')),
          nav_panel('Mean', DTOutput('out_Mu')),
          nav_panel('Covariance Matrix', DTOutput('out_Sigma')),
        )
      )
    )
  )
)

server <- function(input, output, session) {
  output$constrain_items <- renderUI({
    M <- ncol(req(data()))
    K <- tryCatch(round(as.numeric(req(isolate(input$nfactor)))), error = function(e) 0)
    constrains <- paste0('constrain_', 1:M)
    lapply(1:M, function(m) {
      x <- selectInput(constrains[m], paste0('Dimension ', m, ':'), as.list(setNames(1:M, colnames(data()))), m)
      observeEvent(input[[constrains[m]]], {
        select <- sapply(constrains, function(x) as.integer(isolate(input[[x]])))
        k <- which(select == select[m] & (1:M) != m)
        updateSelectInput(session, constrains[k], selected = setdiff(1:M, select))
      })
      if (m <= K)
        x
      else
        hidden(x)
    })
  })
  observeEvent(input$constrain_type, {
    M <- ncol(req(data()))
    label <- paste0(if (input$constrain_type == 'C1') 'Dimension' else 'Dimensions 1 to', ' ', 1:M, ':')
    label[1] <- 'Dimension 1:'
    for (m in 1:M)
      updateSelectInput(session, paste0('constrain_', m), label[m])
  })
  observeEvent(input$nfactor, {
    M <- ncol(req(data()))
    K <- round(as.numeric(req(input$nfactor)))
    if (1 <= K && K <= M) {
      for (m in 1:K)
        showElement(paste0('constrain_', m))
      if (K < M)
        for (m in (K + 1):M)
          hideElement(paste0('constrain_', m))
    }
  })

  data <- reactive({
    dat <- fread(req(input$in_data)$datapath, data.table = F)
    dat[!(dat == 0 | dat == 1)] <- NA
    colnames(dat) <- paste0('Item', 1:ncol(dat))
    updateNumericInput(session, 'nfactor', max = ncol(dat))
    dat
  })
  model <- reactive({
    if (as.integer(input$navbar) < 10)
      return(NULL)
    dat <- fread(req(input$in_model)$datapath, data.table = F)
    dat[is.na(dat)] <- 0
    dat[dat != 0] <- 1
    rownames(dat) <- paste0('Item', 1:nrow(dat))
    colnames(dat) <- paste0('Dim', 1:ncol(dat))
    dat
  })
  group <- reactive({
    if (as.integer(input$navbar) < 20)
      return(NULL)
    raw <- fread(req(input$in_group)$datapath, data.table = F)[, 1]
    raw[is.na(raw)] <- 0
    new <- if (is.integer(raw))
      raw - min(raw) + 1
    else
      as.integer(as.factor(raw))
    data.frame(Raw = raw, New = new)
  })
  output$out_data <- renderDT(formatStyle(datatable(data()), 0, fontWeight = 'bold'))
  output$out_model <- renderDT({
    x <- model()
    if (is.null(x))
      NULL
    else {
      x[x == 0] <- NA
      x[!is.na(x)] <- '✔'
      formatStyle(datatable(x), 0, fontWeight = 'bold')
    }
  })
  output$out_group <- renderDT(group())

  out <- reactiveValues()
  observeEvent(input$navbar, {
    select <- switch(input$navbar,
                     '1' = list(c('GVEM with Post-Hoc Rotation' = 'rotation', 'Lasso GVEM' = 'GVEM', 'Lasso IW-GVEM' = 'IWGVEM', 'Adaptive Lasso GVEM' = 'Adapt GVEM', 'Adaptive Lasso IW-GVEM' = 'Adapt IWGVEM'), 'IWGVEM'),
                     '2' = list(c('SGVEM with Post-Hoc Rotation' = 'rotation', 'Lasso SGVEM' = 'SGVEM', 'Adaptive Lasso SGVEM' = 'Adapt SGVEM'), 'SGVEM'),
                     '11' = list(c('GVEM', 'GVEM-BS', 'GVEM-IW', 'IW-GVEM'), 'IW-GVEM'),
                     '12' = list(c('SGVEM'), 'SGVEM'),
                     '21' = list(c('Lasso EM' = 'EM', 'Lasso EMM' = 'EMM', 'Lasso GVEM' = 'GVEM', 'Lasso IW-GVEMM' = 'IWGVEMM'), 'IWGVEMM'))
    updateSelectInput(inputId = 'algorithm', choices = select[[1]], selected = select[[2]])
    accordion_panel_open('accord_sidebar', c('Data', 'Model'))
  })

  summarize <- function(result, ...) {
    require(tibble)
    a <- as.data.frame(result$a)
    a[a == 0] <- NA
    rownames(a) <- paste0('Item', 1:nrow(a))
    colnames(a) <- paste0('Dim', 1:ncol(a))
    gamma <- result$gamma
    if (!is.null(gamma)) {
      gamma[gamma == 0] <- NA
      d <- dim(gamma)
      if (d[1] == 1)
        gamma <- NULL
      else {
        gamma <- as.data.frame(do.call(rbind, lapply(1:d[1], function(g) {
          x <- gamma[g, , ]
          rownames(x) <- paste0('Group', g, ': Item', 1:nrow(x))
          x
        })))
        colnames(gamma) <- paste0('Dim', 1:ncol(gamma))
      }
    }
    b <- as.data.frame(result$b)
    rownames(b) <- paste0('Item', 1:nrow(b))
    colnames(b) <- paste0('Negative Intercept')
    beta <- result$beta
    if (!is.null(beta)) {
      beta[beta == 0] <- NA
      beta <- as.data.frame(t(beta))
      rownames(beta) <- paste0('Item', 1:nrow(beta))
      colnames(beta) <- paste0('Group', 1:ncol(beta))
    }
    c <- result$c
    if (!is.null(c)) {
      c <- as.data.frame(c)
      rownames(c) <- paste0('Item', 1:nrow(c))
      colnames(c) <- paste0('Guessing')
    }
    Sigma <- result$Sigma
    d <- dim(Sigma)
    if (length(d) == 3 && d[1] > 1)
      Sigma <- do.call(rbind, lapply(1:d[1], function(g) {
        x <- Sigma[g, , ]
        rownames(x) <- paste0('Group', g, ': Dim', 1:nrow(x))
        x
      }))
    else {
      if (length(d) == 3)
        Sigma <- Sigma[1, , ]
      rownames(Sigma) <- paste0('Dim', 1:nrow(Sigma))
    }
    Sigma <- as.data.frame(Sigma)
    colnames(Sigma) <- paste0('Dim', 1:ncol(Sigma))
    Mu <- as.data.frame(if (is.null(result$Mu)) matrix(0, nrow = 1, ncol = ncol(Sigma)) else result$Mu)
    rownames(Mu) <- paste0('Group', 1:nrow(Mu))
    colnames(Mu) <- paste0('Dim', 1:ncol(Mu))
    lst(lambda0 = result$lambda0, a, gamma, b, beta, c, Mu, Sigma, ...)
  }
  estimate <- function(type, method, data, model, group, params) {
    require(VEMIRT)
    require(tibble)
    result <- if (type == '1') {
      if (method == 'rotation') {
        fit <- E2PL_gvem_rot(data, params$nfactor, rot = params$rotation)
        with(fit, list(a = ra, b = rb, Sigma = rsigma))
      } else {
        fit <- if (method %in% c('GVEM', 'IWGVEM'))
          E2PL_gvem_lasso(data, model, constrain = params$constrain, non_pen = params$non_pen)
        else
          E2PL_gvem_adaptlasso(data, model, constrain = params$constrain, non_pen = params$non_pen)
        if (method %in% c('Adapt GVEM', 'Adapt IWGVEM'))
          with(fit, list(a = ra, b = rb, Sigma = rsigma))
        else
          with(E2PL_iw(data, fit), list(a = new_a, b = new_b, Sigma = new_Sigma_theta))
      }
    } else if (type == '2') {
      if (method == 'rotation') {
        fit <- E3PL_sgvem_rot(data, params$nfactor, mu_b = params$b_mu, sigma2_b = params$b_sigma2, Alpha = params$c_alpha, Beta = params$c_beta, rot = params$rotation)
        with(fit, list(a = ra, b = rb, c = rc, Sigma = rsigma))
      } else {
        fit <- if (method == 'SGVEM')
          E3PL_sgvem_lasso(data, model, mu_b = params$b_mu, sigma2_b = params$b_sigma2, Alpha = params$c_alpha, Beta = params$c_beta, constrain = params$constrain, non_pen = params$non_pen)
        else
          E3PL_sgvem_adaptlasso(data, model, mu_b = params$b_mu, sigma2_b = params$b_sigma2, Alpha = params$c_alpha, Beta = params$c_beta, constrain = params$constrain, non_pen = params$non_pen)
        with(fit, list(a = ra, b = rb, c = rc, Sigma = rsigma))
      }
    } else if (type == '11') {
      if (method == 'IW-GVEM') {
        fit <- C2PL_iw2(data, model)
        with(fit$fit, lst(a, b, Sigma))
      } else {
        fit <- C2PL_gvem(data, model)
        if (method == 'GVEM-BS')
          with(C2PL_bs(fit), list(a = boots_a, b = boots_b, Sigma = fit$rsigma))
        else if (method == 'GVEM-IW')
          with(C2PL_iw(data, fit), list(a = new_a, b = new_b, Sigma = new_Sigma_theta))
        else   # GVEM
          with(fit, list(a = ra, b = rb, Sigma = rsigma))
      }
    } else if (type == '12') {
      fit <- C3PL_sgvem(data, model, mu_b = params$b_mu, sigma2_b = params$b_sigma2, Alpha = params$c_alpha, Beta = params$c_beta)
      with(fit, list(a = ra, b = rb, c = rc, Sigma = rsigma))
    } else if (type == '21') {
      Lambda0 <- seq(params$lambda_from, params$lambda_to, by = params$lambda_by)
      fit <- if (method == 'EM' || method == 'EMM')
        D2PL_em(data, model, group, method, Lambda0, params$level)
      else   # Lasso GVEM & Lasso IW-GVEMM
        D2PL_gvem(data, model, group, method, Lambda0)
      with(fit$fit, lst(a, b, Mu = Mu, Sigma = Sigma, gamma = gamma, beta = beta, lambda0 = lambda0))
    } else
      list()
    summarize(result, type, method, data, model, group, params, fit)
  }
  formatDT <- function(x) {
    if (is.null(x))
      NULL
    else {
      row <- nrow(x) > 1
      y <- formatRound(datatable(x, rownames = row), 1:ncol(x), digits = 3, mark = '')
      if (row)
        y <- formatStyle(y, 0, fontWeight = 'bold')
      y
    }
  }
  output$out_a <- renderDT(formatDT(out$a))
  output$out_gamma <- renderDT(formatDT(out$gamma))
  output$out_b <- renderDT(formatDT(out$b))
  output$out_beta <- renderDT(formatDT(out$beta))
  output$out_c <- renderDT(formatDT(out$c))
  output$out_Mu <- renderDT(formatDT(out$Mu))
  output$out_Sigma <- renderDT(formatDT(out$Sigma))
  output$out_lambda0 <- renderText(out$lambda0)

  run <- reactiveValues(label = 'Run', process = r_bg(vector))
  observeEvent(run$label, {
    updateActionButton(session, 'in_run', run$label)
  })
  observeEvent(input$in_run, {
    if (run$label == 'Run') {
      accordion_panel_close('accord_sidebar', 'Data')
      if (as.integer(input$navbar) == 21 && input$algorithm %in% c('EM', 'EMM')) {
        level <- min(max(round(as.numeric(input$level)), 1), 25)
        updateNumericInput(session, 'level', value = level)
      }
      if (as.integer(input$navbar) < 10) {
        nfactor <- min(max(round(as.numeric(input$nfactor)), 1), ncol(isolate(data())))
        updateNumericInput(session, 'nfactor', value = nfactor)
      }
      params <- list(b_mu = as.numeric(input$b_mu), b_sigma2 = as.numeric(input$b_sigma2), c_alpha = as.numeric(input$c_alpha), c_beta = as.numeric(input$c_beta),
                     lambda_from = as.numeric(input$lambda_from), lambda_to = as.numeric(input$lambda_to), lambda_by = as.numeric(input$lambda_by),
                     level = as.integer(input$level), nfactor = as.integer(input$nfactor), rotation = input$rotation, constrain = input$constrain_type)
      args <- lst(type = as.integer(input$navbar), method = input$algorithm, data = data(), model = model(), group = group()$New, params)
      if (args$type < 10 && input$algorithm != 'rotation') {
        M <- ncol(args$data)
        K <- params$nfactor
        id <- diag(rep(1, K))
        if (params$constrain == 'C2')
          id[lower.tri(id)] <- 1
        model <- matrix(1, M, K)
        model[sapply(paste0('constrain_', 1:K), function(x) as.numeric(input[[x]])), ] <- id
        args$model <- model
        args$params$non_pen <- as.numeric(input[[paste0('constrain_', K)]])
      }
      run$process <- r_bg(estimate, args, package = T)
      run$label <- 'Stop'
    } else {
      run$process$kill()
      run$label <- 'Run'
    }
  })
  observe({
    if (run$label != 'Run' && !run$process$is_alive()) {
      hideElement('panel_ic')
      o <- run$process$get_result()
      for (v in names(o))
        out[[v]] <- o[[v]]
      if (out$type == '21') {
        if (out$method %in% c('EM', 'EMM'))
          updateSelectizeInput(session, 'ic', selected = 'BIC')
        else {
          updateSelectizeInput(session, 'ic', selected = 'GIC')
          updateNumericInput(session, 'gic_c', value = 1)
        }
        showElement('panel_ic')
      }
      accordion_panel_set('accord_output', 'Parameter Estimates')
      accordion_panel_open('accord_sidebar', 'Postestimation')
      run$label <- 'Run'
    } else {
      invalidateLater(250)
    }
  })

  observeEvent(list(input$ic, input$gic_c), {
    runjs('Shiny.setInputValue("panel_ic_hidden", $("#panel_ic").is(":hidden"));')
    if (is.null(input$panel_ic_hidden) || input$panel_ic_hidden) return()
    criterion <- if (input$ic == 'GIC')
      as.numeric(input$gic_c)
    else
      input$ic
    o <- summarize(coef(out$fit, criterion))
    for (v in names(o))
      out[[v]] <- o[[v]]
  })
  output$download <- downloadHandler('estimate.xlsx', function(file) {
    add <- function(var, table) {
      if (!is.null(table)) {
        table[is.na(table)] <- 0
        result[[var]] <<- table
      }
    }
    result <- list()
    add('Slope', out$a)
    add('DIF in Slope', out$gamma)
    add('Intercept', out$b)
    add('DIF in Intercept', out$beta)
    add('Guessing', out$c)
    add('Mean', out$Mu)
    add('Covariance Matrix', out$Sigma)
    write.xlsx(result, file, rowNames = T)
  })
}

shinyApp(ui = ui, server = server)
