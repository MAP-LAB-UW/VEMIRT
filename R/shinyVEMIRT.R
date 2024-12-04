# Written by Weicong Lyu

#' Shiny App for VEMIRT
#'
#' @author Weicong Lyu <weiconglyu@um.edu.mo>
#' @export
shinyVEMIRT <- function() {
  suppressWarnings(shiny::runApp(system.file('shiny', package = 'VEMIRT'), display.mode = 'normal'))
}
