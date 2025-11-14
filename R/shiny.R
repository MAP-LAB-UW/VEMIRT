#' Shiny App for VEMIRT
#'
#' @author Weicong Lyu <weiconglyu@um.edu.mo>
#' @export
shinyVEMIRT <- function() {
  suppressWarnings(shiny::runApp(system.file('shinyVEMIRT', package = 'VEMIRT'), display.mode = 'normal'))
}

#' Shiny App for DIF Dashboard
#'
#' @author Yijun Cheng <chengxb@uw.edu>
#' @author He Ren <heren@uw.edu>
#' @author Weicong Lyu <weiconglyu@um.edu.mo>
#' @export
DIFdashboard <- function() {
  suppressWarnings(shiny::runApp(system.file('DIFdashboard', package = 'VEMIRT'), display.mode = 'normal'))
}
