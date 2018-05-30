#' Georacle
#'
#' When the function is invoked the Goeracle shiny app is started. 
#' 
#' @return
#' Does not return anything
#' @examples
#' # The following commands will load the shiny app either through an RStudio session or
#' # through your internet browser
#'
#' library("Georacle")
#' \donttest{ Georacle() }
#'
#'
#' @export
Georacle <- function()
{
  appDir <- system.file("shiny-app",  package = "Georacle")

  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `Georacle`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode="normal")

}
