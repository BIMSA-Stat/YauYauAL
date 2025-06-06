#' Run the YYF Shiny Application
#'
#' A convenience function to launch the Shiny app from this package.
#' @export
run_app <- function() {
  appDir <- system.file("shinyApp", package = "YauYauAL")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `YauYauAL`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
