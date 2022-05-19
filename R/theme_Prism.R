#' Theme for imitating Prism
#'
#' @param base_size Numerical. Font size.
#' @param base_family Character. Font family.
#'
#' @return ggplot2 theme
#'
#' @export


theme_Prism <- function (
  base_size=7,
  base_family="Arial"
  ){

  half_line <- base_size/2
  ggplot2::theme(
    #1.436 in illustrator is 1.436/(ggplot2::.pt*72.27/96); 0.67 is a close approximation
    line = ggplot2::element_line(colour = "black", size = 0.67, linetype = 1, lineend = "square"),
    rect = ggplot2::element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1),
    text = ggplot2::element_text(family = base_family, face = "bold", colour = "black", size = base_size,
                        lineheight = 0.9,  hjust = 0.5, vjust = 0.5, angle = 0, margin = ggplot2::margin(), debug = FALSE),

    axis.line = ggplot2::element_line(colour = "black", size = 0.67, linetype = 1, lineend = "square"),
    axis.line.y.left = ggplot2::element_line(colour = "black", size = 0.67, linetype = 1, lineend = "square"),
    axis.text = ggplot2::element_text(size = base_size, colour = "black", face = "bold", family = base_family),
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0.8*half_line/2), face = "bold", family = base_family,
                               vjust = 1),
    axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 0.8*half_line/2), face = "bold", family = base_family,
                               hjust = 1),
    axis.ticks = ggplot2::element_line(colour = "grey20"),
    #axis.ticks.length = unit(half_line/2, "pt"),
    axis.ticks.length = grid::unit(3.9, "pt"),
    axis.title.x = ggplot2::element_text(face = "bold", family = base_family,
                                margin = ggplot2::margin(t = 0.8 * half_line,
                                                b = 0.8 * half_line/2)),
    axis.title.y = ggplot2::element_text(angle = 90, face = "bold", family = base_family,
                                margin = ggplot2::margin(r = 0.8 * half_line,
                                                l = 0.8 * half_line/2)),

    legend.background = ggplot2::element_blank(),
    legend.margin = ggplot2::margin(t=0.2, unit="cm"),
    legend.key = ggplot2::element_blank(),
    legend.key.size = grid::unit(1.2, "lines"),
    legend.key.height = NULL,
    legend.key.width = NULL,
    legend.text = ggplot2::element_text(size = base_size, face = "bold", family = base_family),
    legend.text.align = NULL,
    legend.title = ggplot2::element_text(hjust = 0, face = "bold", family = base_family),
    legend.title.align = NULL,
    legend.position = "right",
    legend.direction = NULL,
    legend.justification = "center",
    legend.box = NULL,

    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing = grid::unit(7, "pt"),
    panel.margin.x = NULL,
    panel.spacing.y = NULL,
    panel.ontop = TRUE,

    strip.background = ggplot2::element_rect(fill = "white", colour = "black", size = 1, linetype = 1),
    strip.text = ggplot2::element_text(colour = "black", size = base_size, face = "bold", family = base_family),
    strip.text.x = ggplot2::element_text(face = "bold", family = base_family,
                                margin = ggplot2::margin(t = half_line,
                                                b = half_line)),
    strip.text.y = ggplot2::element_text(angle = -90, face = "bold", family = base_family,
                                margin = ggplot2::margin(l = half_line,
                                                r = half_line)),
    strip.switch.pad.grid = grid::unit(0.1, "cm"),
    strip.switch.pad.wrap = grid::unit(0.1, "cm"),

    plot.background = ggplot2::element_rect(colour = "white"),
    plot.title = ggplot2::element_text(size = base_size, face = "bold", family = base_family,
                              margin = ggplot2::margin(b = half_line * 1.2)),
    plot.margin = ggplot2::margin(0.1,0.1,0.1,0.1),
    complete = TRUE
  )
}
