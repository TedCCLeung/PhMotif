#' Function to plot motif enrichment results as heat maps
#'
#' @param dir Character. Output directory.
#' @param enrichment_dir Character. Full path to the run_enrich output directory.
#' @param motif_MEME_file Character. Full path to a motif .meme file.
#' @param enrich_pval_threshold Numerical. Default 1e-2. P value threshold to decide for an enriched motif.
#' @param cutree_height Numerical. Similarity value to cut the dendrogram for merging motifs.
#' @param remove_zoom whether to remove zoom box
#'
#' @return None.

run_heatmap <- function(
  dir,
  enrichment_dir,
  motif_MEME_file,
  enrich_pval_threshold = 1e-2,
  cutree_height = 0.5,
  remove_zoom = FALSE
){

  ## STEP 0: SUMMARISE RESULTS -----------------------------
  if (!dir.exists(paste0(dir, "/summary/"))){dir.create(paste0(dir, "/summary/"), recursive = TRUE)}
  if (!dir.exists(paste0(dir, "/heatmaps/"))){dir.create(paste0(dir, "/heatmaps/"), recursive = TRUE)}

  ## Retrieve the data from fisher enrichment test ----
  enrich_result_summary(dir = enrichment_dir, max_pval = enrich_pval_threshold, outdir = paste0(dir, "/summary/"))

  ## Read in information ---
  Fisher <- read_Fisher_result(paste0(dir, "/summary/"))
  motifs <- motif_MEME_file %>% universalmotif::read_meme() %>% universalmotif::filter_motifs(name = Fisher$motif_names)

  ## Clustering ---
  hclust_out <- motifs %>% hclust_motifs()
  cutree_out <- stats::cutree(hclust_out, h = cutree_height)

  ## Heat maps ---
  heatmap_PVAL <- heatmap_from_Fisher(Fisher$mat_logPVAL, hclust_out, cutree_out, motifs, "PVAL", remove_zoom = remove_zoom)
  heatmap_FDR <- heatmap_from_Fisher(Fisher$mat_logFDR, hclust_out, cutree_out, motifs, "FDR")
  heatmap_FC <- heatmap_from_Fisher(Fisher$mat_FC, hclust_out, cutree_out, motifs, "FC")

  ## Get the motif families for plotting the bar ----
  ## Any of the heatmap can be used for the main
  heatmap_TF_family <- heatmap_motif_to_CISBP_family(Fisher$mat_FC, hclust_out, motifs, main = heatmap_FC) %>% heatmap_TF_family()
  heatmap_motif_source <- heatmap_motif_source(matrix = Fisher$mat_FC, hclust_out = hclust_out, main = heatmap_FC)

  ## Add the extra information to the heat maps ----
  heatmap_PVAL <- heatmap_PVAL + heatmap_TF_family + heatmap_motif_source
  heatmap_FDR <- heatmap_FDR + heatmap_TF_family + heatmap_motif_source
  heatmap_FC <- heatmap_FC + heatmap_TF_family + heatmap_motif_source

  ## Output pdf ----

  grDevices::pdf(paste0(dir, "heatmaps/heatmap_PVAL.pdf"))
  ComplexHeatmap::draw(heatmap_PVAL)
  grDevices::dev.off()

  grDevices::pdf(paste0(dir, "heatmaps/heatmap_FDR.pdf"))
  ComplexHeatmap::draw(heatmap_FDR)
  grDevices::dev.off()

  grDevices::pdf(paste0(dir, "heatmaps/heatmap_FC.pdf"))
  ComplexHeatmap::draw(heatmap_FC)
  grDevices::dev.off()
}
