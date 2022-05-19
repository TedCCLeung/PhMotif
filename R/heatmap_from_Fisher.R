#' Plot results of Fisher Exact Test as heat maps
#'
#' @param matrix Matrix. Numerical matrix to be plotted
#' @param hclust_out hclust object.
#' @param cutree_out cutree object.
#' @param motifs Universalmotifs object.
#' @param style Character. One of "FC", "PVAL" and "FDR".
#' @param remove_zoom whether to remove zoom box
#'
#' @return ComplexHeatmap object
#'
#' @export

heatmap_from_Fisher <- function(
  matrix,
  hclust_out,
  cutree_out,
  motifs,
  style,
  remove_zoom = FALSE
){

  if (!style %in% c("FC", "PVAL", "FDR")){stop("Wrong style entered!")}

  if (style == "FC"){heatmap <- heatmap_FC(matrix, hclust_out, cutree_out, motifs, remove_zoom = remove_zoom)
  } else if (style == "PVAL"){heatmap <- heatmap_PVAL(matrix, hclust_out, cutree_out, motifs, remove_zoom = remove_zoom)
  } else if (style == "FDR"){heatmap <- heatmap_FDR(matrix, hclust_out, cutree_out, motifs, remove_zoom = remove_zoom)
  }

  return(heatmap)
}



heatmap_FC <- function(
  matrix,
  hclust_out,
  cutree_out,
  motifs,
  remove_zoom = FALSE
){

  heatmap <- suppressMessages(
    heatmap_motif_clustering(
      input_matrix = matrix,
      hclust_out = hclust_out,
      cutree_out = cutree_out,
      motifset = motifs,
      colors = c("#3B4992", "#FFFFFF", "#A20056"),
      color_range = c(-2, 0, 2),
      scale_name = " ",
      panel_size = grid::unit(0.3, "cm"),
      panel_gap = grid::unit(0.0, "cm"),
      panel_width = grid::unit(2.0, "cm"),
      panel_link_width = grid::unit(0.6, "cm"),
      remove_zoom = remove_zoom
      )
    )

  return(heatmap)
}


heatmap_PVAL <- function(
  matrix,
  hclust_out,
  cutree_out,
  motifs,
  remove_zoom = FALSE
){

  heatmap <- suppressMessages(
    heatmap_motif_clustering(
      input_matrix = matrix,
      hclust_out = hclust_out,
      cutree_out = cutree_out,
      motifset = motifs,
      colors = c("#EE0000", "#FFFFFF"),
      color_range = c(6, 0),
      scale_name = " ",
      panel_size = grid::unit(0.3, "cm"),
      panel_gap = grid::unit(0.0, "cm"),
      panel_width = grid::unit(2.0, "cm"),
      panel_link_width = grid::unit(0.6, "cm"),
      remove_zoom = remove_zoom
    )
  )

  return(heatmap)
}


heatmap_FDR <- function(
  matrix,
  hclust_out,
  cutree_out,
  motifs,
  remove_zoom = FALSE
){

  heatmap <- suppressMessages(
    heatmap_motif_clustering(
      input_matrix = matrix,
      hclust_out = hclust_out,
      cutree_out = cutree_out,
      motifset = motifs,
      colors = c("#EE0000", "#FFFFFF"),
      color_range = c(6, 0),
      scale_name = " ",
      panel_size = grid::unit(0.3, "cm"),
      panel_gap = grid::unit(0.0, "cm"),
      panel_width = grid::unit(2.0, "cm"),
      panel_link_width = grid::unit(0.6, "cm"),
      remove_zoom = remove_zoom
    )
  )

  return(heatmap)
}




