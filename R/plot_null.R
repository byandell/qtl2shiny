#' @importFrom ggplot2 ggplot aes geom_text theme_void
plot_null <- function() {
  ggplot2::ggplot(data.frame(x=1,y=1), 
                  ggplot2::aes(x,y,label="no data")) +
    ggplot2::geom_text(size=10) + 
    ggplot2::theme_void()
}
