#' Function to take in motifs and sequences as inputs
#'
#' @param dir Character. Output directory.
#' @param motifs Either a character vector of full paths to .meme motifs files, or a list of universalmotif objects.
#' @param sequences Either a character vector of full paths to .fa DNA sequence files, of a list of Biostring DNAStringSet / DNAString objects.
#'
#' @return None.


run_input <- function(
  dir,
  motifs,
  sequences
){

  ## STEP 0: SET UP DIRECTORIES -----------------------------

  if (!dir.exists(paste0(dir, "/MEME"))){dir.create(paste0(dir, "/MEME"), recursive = TRUE)}
  if (!dir.exists(paste0(dir, "/FASTA"))){dir.create(paste0(dir, "/FASTA"), recursive = TRUE)}

  ## STEP 1: Process motifs -----------------------------

  motifs <- input_motifs(motifs, out_file = paste0(dir, "/MEME/motifs.meme"))

  ## STEP 2: Process sequences -----------------------------

  if (class(unlist(sequences[[1]])) == "character"){
    sequences <- sequences %>% lapply(read_genes) %>%
      genes_to_promoter_Biostring(rDEI_threshold = 1, upstream = 1000, downstream = 500)
    file_names <- get_file_names(sequences)
  } else if (class(unlist(sequences)[[1]]) %in% c("DNAStringSet", "DNAString")){
    sequences <- sequences
    file_names <- paste0(names(sequences), ".fa")
  } else {
    stop("Input the correct sequences! Either Biostrings or path to files!")
  }

  for (k in 1:length(file_names)){Biostrings::writeXStringSet(sequences[[k]], filepath = paste0(dir, "FASTA/", file_names[[k]]))}
}


