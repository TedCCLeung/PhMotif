#' Function to perform motif mapping
#'
#' @param dir Character. Output directory.
#' @param FASTA_files Character vector of full paths to FASTA .fa files.
#' @param motif_MEME_file Character. Full path to a motif .meme file.
#' @param mapping_pval_threshold Numerical. Default 1e-4. Threshold to determine whether a motif has been mapped to a stretch of sequence.
#'
#' @return None.

run_mapping <- function(
  dir,
  FASTA_files,
  motif_MEME_file,
  mapping_pval_threshold = 1e-4
){

  ## STEP 0: SET UP DIRECTORY -----------------------------
  if (!dir.exists(dir)){dir.create(dir, recursive = TRUE)}

  ## STEP 1: SET UP FUNCTION FOR INDIVIDUAL MAPPINGS -----------------------------
  motifs <- universalmotif::read_meme(motif_MEME_file)

  mapping_motifs <- function(
    FASTA_file,
    file_name
  ){
    scan_seq_out <- universalmotif::scan_sequences(
      motifs = motifs,
      sequences = Biostrings::readDNAStringSet(FASTA_file),
      RC = TRUE,
      threshold = mapping_pval_threshold,
      threshold.type = "pvalue",
      no.overlaps = TRUE,
      no.overlaps.strat = "score",
      return.granges = TRUE,
      calc.qvals.method = "BH"
    )
    rtracklayer::export(scan_seq_out, file_name, format = "gff3")
  }


  ## STEP 2: LOOP THEM WITH MAPPLY -----------------------------

  sample_names <- FASTA_files %>%
    strsplit(split = "/") %>%
    lapply(utils::tail, 1) %>%
    base::unlist()
  file_names <- paste0(dir, sample_names, ".gff3")

  result <- base::mapply(mapping_motifs, FASTA_files, file_names)
}

