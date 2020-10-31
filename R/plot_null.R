#' @importFrom ggplot2 ggplot aes geom_text theme_void
#' @importFrom rlang .data
#' 
plot_null <- function(msg = "no data") {
  ggplot2::ggplot(data.frame(x = 1, y = 1), 
                  ggplot2::aes(.data$x, .data$y, label = msg)) +
    ggplot2::geom_text(size = 10) + 
    ggplot2::theme_void()
}
