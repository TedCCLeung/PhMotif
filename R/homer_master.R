homer_master <- function(
  motifs = universalmotif::read_homer("/Users/TedCCLeung/Documents/Projects/Packages/all.8.motif"),
  dir = "/Users/TedCCLeung/Documents/Projects/Photoperiod/2_analysis/2_pipeline/PhotoperiodMotif/homer_08/",
  known_only = TRUE,
  pval_threshold = 1e-3,
  FC_threshold = 0.5,
  cutree_height = 0.55
){

  df_summary <- input_homer(dir)
  utils::write.table(df_summary, paste0(dir, "summary_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  if (known_only){
    df_summary <- df_summary[startsWith(df_summary$name, "M"), ]
  }

  heatmap_list <- process_homer_matrix(
    dir = dir,
    df_result = df_summary,
    motifs = motifs,
    pval_threshold = pval_threshold,
    FC_threshold = FC_threshold,
    cutree_height = cutree_height
  )

  grDevices::pdf(paste0(dir, "heatmap_PVAL.pdf"), height = 20, width = 5.5)
  ComplexHeatmap::draw(heatmap_list$heatmap_PVAL)
  grDevices::dev.off()

  grDevices::pdf(paste0(dir, "heatmap_FDR.pdf"), height = 20, width = 5.5)
  ComplexHeatmap::draw(heatmap_list$heatmap_FDR)
  grDevices::dev.off()

  grDevices::pdf(paste0(dir, "heatmap_FC.pdf"), height = 20, width = 5.5)
  ComplexHeatmap::draw(heatmap_list$heatmap_FC)
  grDevices::dev.off()
}


# homer_master(
#   motifs = universalmotif::read_homer("/Users/TedCCLeung/Documents/Projects/Packages/all.8.motif"),
#   dir = "/Users/TedCCLeung/Documents/Projects/Photoperiod/2_analysis/2_pipeline/PhotoperiodMotif/homer_08/",
#   pval_threshold = 5e-3,
#   FC_threshold = 0.5,
#   cutree_height = 0.55
# )
#
# homer_master(
#   motifs = universalmotif::read_homer("/Users/TedCCLeung/Documents/Projects/Packages/all.8.motif"),
#   dir = "/Users/TedCCLeung/Documents/Projects/Photoperiod/2_analysis/2_pipeline/PhotoperiodMotif/homer_10/",
#   pval_threshold = 5e-4,
#   FC_threshold = 1.0,
#   cutree_height = 0.55,
#   known_only = FALSE
# )
