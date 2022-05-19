process_homer_matrix <- function(
  dir,
  df_result,
  motifs,
  pval_threshold = 1e-4,
  FC_threshold = 1, # log
  cutree_height = 0.5
){

  ## Make the matrices
  mat_pval <- df_result[, c("name", "sample", "pvalue")] %>%
    tidyr::pivot_wider(names_from = "sample", values_from = "pvalue") %>%
    tibble::column_to_rownames(var = "name") %>%
    log10() %>%
    as.matrix() %>%
    dplyr::na_if(Inf) %>%
    dplyr::na_if(-Inf) %>%
    tidyr::replace_na(0)

  mat_fdr <- df_result[, c("name", "sample", "fdr")] %>%
    tidyr::pivot_wider(names_from = "sample", values_from = "fdr") %>%
    tibble::column_to_rownames(var = "name") %>%
    log10() %>%
    as.matrix() %>%
    dplyr::na_if(Inf) %>%
    dplyr::na_if(-Inf) %>%
    tidyr::replace_na(0)

  mat_fc <- df_result[, c("name", "sample", "fold_change")] %>%
    tidyr::pivot_wider(names_from = "sample", values_from = "fold_change") %>%
    tibble::column_to_rownames(var = "name") %>%
    suppressWarnings() %>%
    as.matrix() %>%
    dplyr::na_if(Inf) %>%
    dplyr::na_if(-Inf) %>%
    tidyr::replace_na(0)

  ## Filter
  row_min <- apply(mat_pval, 1, min)
  FC_max <- apply(mat_fc, 1, max)
  common_motifs <- rownames(mat_pval)[(row_min < log10(pval_threshold)) & (FC_max > FC_threshold)]
  motifs <- universalmotif::filter_motifs(motifs, name = common_motifs)

  mat_pval <- mat_pval[rownames(mat_pval) %in% common_motifs, ]*-1
  mat_fdr <- mat_fdr[rownames(mat_fdr) %in% common_motifs, ]*-1
  mat_fc <- mat_fc[rownames(mat_fc) %in% common_motifs, ]

  ## Clustering ---
  hclust_out <- motifs %>% hclust_motifs() %>% suppressWarnings()
  cutree_out <- stats::cutree(hclust_out, h = cutree_height)

  ## Heat maps ---
  heatmap_PVAL <- heatmap_from_Fisher(mat_pval, hclust_out, cutree_out, motifs, "PVAL")
  heatmap_FDR <- heatmap_from_Fisher(mat_fdr, hclust_out, cutree_out, motifs, "FDR")
  heatmap_FC <- heatmap_from_Fisher(mat_fc, hclust_out, cutree_out, motifs, "FC")

  ## Get the motif families for plotting the bar ----
  ## Any of the heatmap can be used for the main
  #heatmap_TF_family <- heatmap_motif_to_CISBP_family(matrix = mat_fc, hclust_out, motifs, main = heatmap_FC) %>% heatmap_TF_family() %>% suppressWarnings()
  heatmap_motif_source <- heatmap_motif_source(matrix = mat_fc, hclust_out = hclust_out, main = heatmap_FC) %>% suppressWarnings()

  ## Add the extra information to the heat maps ----
  heatmap_PVAL <- heatmap_PVAL + heatmap_motif_source #+ heatmap_TF_family
  heatmap_FDR <- heatmap_FDR + heatmap_motif_source #+ heatmap_TF_family
  heatmap_FC <- heatmap_FC + heatmap_motif_source #+ heatmap_TF_family

  ## Output pdf ----
  return(list(heatmap_PVAL = heatmap_PVAL, heatmap_FDR = heatmap_FDR, heatmap_FC = heatmap_FC))
}
