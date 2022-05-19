#' Function to perform motif enrichment pipeline from motifs and sequences.
#'
#' @param output_dir Character. Output directory.
#' @param motifs Either a character vector of full paths to .meme motifs files, or a list of universalmotif objects.
#' @param sequences Either a character vector of full paths to .fa DNA sequence files, of a list of Biostring DNAStringSet / DNAString objects.
#' @param skip_input skip input
#' @param skip_mapping skip mapping
#' @param skip_enrich skip enrich
#' @param skip_heatmap skip heat map
#'
#' @return None.


run_master <- function(
  output_dir,
  motifs,
  sequences,
  skip_input = FALSE,
  skip_mapping = FALSE,
  skip_enrich = FALSE,
  skip_heatmap = FALSE
){

  ## STEP 0: PROCESS INPUT MOTIFS AND SEQUENCES -----------------------------
  if (!skip_input){
    run_input(
      dir = paste0(output_dir, "/0_input/"),
      motifs = motifs,
      sequences = sequences
    )
  }

  ## STEP 1: MAPPING OF THE MOTIFS ONTO PROMOTERS -----------------------------
  if (!skip_mapping){
    run_mapping(
      dir = paste0(output_dir, "/1_mapping/"),
      FASTA_files = list.files(paste0(output_dir, "/0_input/FASTA"), full.names = TRUE),
      motif_MEME_file = paste0(output_dir, "/0_input/MEME/motifs.meme"),
      mapping_pval_threshold = 1e-4
    )
  }

  ## STEP 2: FISHER ENRICHMENT TESTS -----------------------------
  if (!skip_enrich){
    run_enrich(
      dir = paste0(output_dir, "/2_enrichment"),
      FASTA_files = list.files(paste0(output_dir, "/0_input/FASTA/"), full.names = TRUE),
      GFF3_files = list.files(paste0(output_dir, "/1_mapping/"), full.names = TRUE)
    )
  }

  ## STEP 3: PLOT HEATMAPS -----------------------------
  if (!skip_heatmap){
    run_heatmap(
      dir = paste0(output_dir, "/3_heatmap/"),
      motif_MEME_file = paste0(output_dir, "/0_input/MEME/motifs.meme"),
      enrichment_dir = paste0(output_dir, "/2_enrichment/results/"),
      enrich_pval_threshold = 1,
      cutree_height = 0.35,
      remove_zoom = FALSE
    )
  }

}

## Both
# run_master(
#   output_dir = "/Users/TedCCLeung/Documents/Projects/Photoperiod/2_analysis/2_pipeline/PhotoperiodMotif/static/",
#   motifs = universalmotif::filter_motifs(motifs = c(CISBP_motifset, denovo_motifset1), name = get_motif_names_from_RF()),
#   sequences = example_promoters
# )



