heatmap_motif_source <- function(
  matrix,
  hclust_out,
  main
){

  ## GET LABELS --------------------------------

  ## Get the motif order in the main heatmap
  #main_row_lab <- rownames(matrix[hclust_out$labels,])[unlist(ComplexHeatmap::row_order(main))]
  main_row_lab <- rownames(matrix[hclust_out$labels,])

  CISBP_motif_names <- CISBP_families$motif

  name_df <- universalmotif::filter_motifs(CISBP_motifset, name = main_row_lab) %>% universalmotif::to_df()
  name_vec <- name_df$altname
  names(name_vec) <- name_df$name

  main_row_lab_ <- main_row_lab
  main_row_lab_[main_row_lab_ %in% CISBP_motif_names] <- paste0(main_row_lab_[main_row_lab_ %in% CISBP_motif_names], "-", name_vec[main_row_lab_[main_row_lab_ %in% CISBP_motif_names]])

  source <- rep("CISBP", length(main_row_lab))
  source[startsWith(main_row_lab, "weeder")] <- "weeder"
  source[startsWith(main_row_lab, "streme")] <- "streme"
  source[startsWith(main_row_lab, "homer")] <- "homer"
  names(source) <- main_row_lab

  mat_source <- as.matrix(source)

  bar <- ComplexHeatmap::Heatmap(
    ## Input
    mat_source,
    show_row_names = TRUE,
    cluster_rows = FALSE,
    row_labels = main_row_lab_,
    ## Color
    col = c("streme" = "#A20056", "CISBP" = "#595959", "weeder" = "#008B45", "homer" = "#008280"),
    name = "Source"
  )

  ## RETURN RESULTS --------------------------------
  return(bar)
}
