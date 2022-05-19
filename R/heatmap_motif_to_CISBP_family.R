#' Retrieve the motifs from a motif heatmap
#'
#' @param matrix input_matrix of the main heatmap
#' @param hclust_out hclust_out used to generate of the main heatmap
#' @param main ComplexHeatmap object.
#' @param motifs Universalmotif object. The same one used to generate the main heatmap
#' @param motif_similarity_cutoff similarity cutoff to determine closest motif in CISBP
#'
#' @return Character vector
#'
#' @export

heatmap_motif_to_CISBP_family <- function(
  matrix,
  hclust_out,
  motifs,
  main,
  motif_similarity_cutoff = 0.9
){

  ## TF FAMILY --------------------------------

  ## Get the motif order in the main heatmap
  main_row_lab <- rownames(matrix[hclust_out$labels,])[unlist(ComplexHeatmap::row_order(main))]

  ## Make motif df
  final_motifs <- universalmotif::filter_motifs(motifs = motifs, name = main_row_lab) %>% universalmotif::to_df()
  ## Make the motifs in order then turn back into list
  final_motifs <- final_motifs[match(main_row_lab, final_motifs$name), ] %>% universalmotif::to_list()

  motif_family <- motif_to_CISBP_family(final_motifs, similarity_cutoff = motif_similarity_cutoff)

  return(motif_family)
}
