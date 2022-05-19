#' Plot the result of hierarchical clustering of motifs
#'
#' @param input_matrix Numerical matrix with motif names as row names and sample names as column names.
#' @param motifset Universalmotif object
#' @param hclust_out result from hclust() of the motifs
#' @param cutree_out result from cutree()
#' @param scale_name Title of the gradient scale legend.
#' @param colors Character vector of colors for the gradient.
#' @param color_range Numerical vector that represents the values the number correspond to. Must be of the same length as colors.
#' @param panel_size Size of the panel in ggplot::grid::unit().
#' @param panel_gap Gap of the panel in ggplot::grid::unit().
#' @param panel_width Width of the panel in ggplot::grid::unit().
#' @param panel_link_width Link width of the panel in ggplot::grid::unit().
#'
#' @return ComplexHeatmap object
#'
#' @export

heatmap_motif_clustering <- function(
  input_matrix,
  hclust_out,
  cutree_out,
  motifset,
  colors = c("purple", "white", "orange"),
  color_range = c(-2, 0, 2),
  scale_name = " ",
  panel_size = grid::unit(0.75, "cm"),
  panel_gap = grid::unit(0.0, "cm"),
  panel_width = grid::unit(4.0, "cm"),
  panel_link_width = grid::unit(0.6, "cm"),
  remove_zoom = FALSE
){

  ## APPERANCE --------------------------------

  ## Options for the heat map appearance
  ComplexHeatmap::ht_opt(
    legend_title_gp = grid::gpar(fontsize = 7, fontface = "bold", base_family="Arial"),
    legend_labels_gp = grid::gpar(fontsize = 7, fontface = "bold", base_family="Arial"),
    heatmap_column_names_gp = grid::gpar(fontsize = 7, fontface = "bold", base_family="Arial"),
    heatmap_column_title_gp = grid::gpar(fontsize = 7, fontface = "bold", base_family="Arial"),
    heatmap_row_title_gp = grid::gpar(fontsize = 7, fontface = "bold", base_family="Arial"),
    heatmap_row_names_gp = grid::gpar(fontsize = 7, fontface = "bold", base_family="Arial"),
    message = FALSE
  )
  color_gradient = circlize::colorRamp2(color_range, colors)

  ## SET PANELS --------------------------------

  ## Panel function
  panel_fun <- function(
    index
  ){
    motif_names <- rownames(input_matrix[hclust_out$labels,])[index]
    merged_motif <- motifset %>%
      universalmotif::filter_motifs(name = motif_names) %>%
      universalmotif::merge_motifs() %>%
      universalmotif::trim_motifs()

    g <- universalmotif::view_motifs(merged_motif) +
      cowplot::theme_nothing()
    g <- grid::grid.grabExpr(print(g))

    grid::pushViewport(grid::viewport(gp = grid::gpar(col = "#FFFFFF")))
    grid::grid.rect()
    grid::grid.draw(g)
    grid::popViewport()
  }

  ## Get the alignment of the panels to the heatmap
  panel_order <- ComplexHeatmap::Heatmap(
    input_matrix[hclust_out$labels,],
    cluster_rows = dendextend::color_branches(hclust_out, k = max(cutree_out)),
    row_split = max(cutree_out)
  ) %>%
    ComplexHeatmap::draw() %>%
    ComplexHeatmap::row_order()

  ## Annotation function
  anno <- ComplexHeatmap::anno_zoom(
    align_to = panel_order,
    panel_fun = panel_fun,

    which = "row",
    size = panel_size,

    gap = panel_gap,
    link_gp = grid::gpar(col = "#CDCDCD"),
    link_width = panel_link_width,

    width = panel_width,



  )

  ## MAIN HEATMAP ACTUALLY OUTPUTTED --------------------------------

  if (remove_zoom){

    main <- ComplexHeatmap::Heatmap(
      ## Input
      input_matrix[hclust_out$labels,],
      ## Columns
      column_order = colnames(input_matrix),
      cluster_column_slices = FALSE,
      ## Rows
      cluster_rows = dendextend::color_branches(hclust_out, k = max(cutree_out)),
      row_split = max(cutree_out),
      #show_row_names = TRUE,
      #right_annotation = ComplexHeatmap::rowAnnotation(foo = anno),
      ## Color
      use_raster = TRUE,
      name = scale_name,
      col = color_gradient
    )

  } else {

    main <- ComplexHeatmap::Heatmap(
      ## Input
      input_matrix[hclust_out$labels,],
      ## Columns
      column_order = colnames(input_matrix),
      cluster_column_slices = FALSE,
      ## Rows
      cluster_rows = dendextend::color_branches(hclust_out, k = max(cutree_out)),
      row_split = max(cutree_out),
      #show_row_names = TRUE,
      right_annotation = ComplexHeatmap::rowAnnotation(foo = anno),
      ## Color
      use_raster = TRUE,
      name = scale_name,
      col = color_gradient
    )

  }



  return(main)
}
