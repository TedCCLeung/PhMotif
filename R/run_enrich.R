#' Function to perform motif enrichment from motif mapping result
#'
#' @param dir Character. Output directory.
#' @param FASTA_files Character vector of full paths to FASTA .fa files.
#' @param GFF3_files Character vector of full paths to GFF3 .gff3 files.
#'
#' @return None.

run_enrich <- function(
  dir,
  FASTA_files,
  GFF3_files
){

  ## STEP 0: SET UP DIRECTORIES -----------------------------
  if (!dir.exists(paste0(dir, "/background"))){dir.create(paste0(dir, "/background"), recursive = TRUE)}
  if (!dir.exists(paste0(dir, "/results"))){dir.create(paste0(dir, "/results"), recursive = TRUE)}


  ## STEP 1: GENERATE BACKGROUNDS -----------------------------

  ## SEQUENCES ---
  concatenate_FASTA(
    input_FASTA =  FASTA_files,
    output_FASTA = paste0(dir, "/background/background.fa")
  )

  ## MAPPINGS ---
  concatenate_gff3(
    input_gff3 = GFF3_files,
    output_gff3 = paste0(dir, "/background/mapping.gff3")
  )


  ## STEP 2: PERFORM MOTIF ENRICHMENT -----------------------------

  ## Set up function to perform individual motif enrichment tests
  fun <- function(
    GFF3_file,
    FASTA_file,
    output_file
  ){
    fisher_result <- enrich_motifs_fisher(
      qry_gff3 = GFF3_file,
      qry_fa = FASTA_file,
      bkg_gff3 = paste0(dir, "/background/mapping.gff3"),
      bkg_fa = paste0(dir, "/background/background.fa"),
      alternative = "two.sided",
      motif_max_p = 1e-5
    )
    utils::write.table(fisher_result, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  ## Loop over the files with mapply
  sample_names <-  FASTA_files %>% strsplit(split = "/") %>% sapply(utils::tail, 1) %>% strsplit(split = "\\.") %>% sapply(utils::head, 1)
  file_names <- paste0(dir, "/results/", sample_names, ".tsv")
  mapply(fun, GFF3_files, FASTA_files, file_names)
}
