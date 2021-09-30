#' Popup shiny function from Joe Cheng
#' 
#' @param title Title of help popup
#' @param content Content of help as character string
#' @param placement Placement position
#' @param trigger Trigger mechanism
#' 
#' @details{
#' See <https://gist.github.com/jcheng5/5913297>.
#' }
#' 
#' @importFrom shiny tagList singleton
#' 
helpPopup <- function(title, content,
                      placement = c('right', 'top', 'left', 'bottom'),
                      trigger = c('click', 'hover', 'focus', 'manual')) {
  shiny::tagList(
    shiny::singleton(
      tags$head(
        tags$script("$(function() { $(\"[data-toggle='popover']\").popover()})"),
        tags$style(type = "text/css",
                   ".popover{max-width:500px; position: fixed; color: black;}")
      )
    ),
    tags$a( # Trying to make title have non-white background with no success.
      href = "#",
      class = "btn btn-link",
      `data-toggle` = "popover", `data-html` = "true",
      title = title,
      `data-content` = content,
      `data-animation` = TRUE,
      `data-placement` = match.arg(placement, several.ok = TRUE)[1],
      `data-trigger` = match.arg(trigger, several.ok = TRUE)[1],
      "Help ..."
    )
  )
}
