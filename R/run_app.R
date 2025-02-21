#' Run the YYF Shiny Application
#'
#' A convenience function to launch the Shiny app from this package.
#' @export
run_app <- function() {
  appDir <- system.file("shinyApp", package = "YauYauFilter")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `YauYauFilter`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
