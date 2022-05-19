#' Plot the result from analyze_positional_bias
#'
#' @param df Output from the analyze_positional_bias function
#'
#' @return ggplot2 object
#'
#' @export

plot_positional_bias <- function(
  df
){
  p <- ggplot2::ggplot(data = df, ggplot2::aes_string(x = "x", y = "y", color = "sequence")) +
    ggplot2::geom_line(ggplot2::aes_string(color = "sequence")) +
    theme_Prism()
  return(p)
  #ggplot2::ggsave(p, filename = outfile, height = 2.5, width = 4)    ## save the figure
}
