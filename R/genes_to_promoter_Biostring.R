#' Function to filter gene IDs (TAIR IDs) based on daily expression integral
#'
#' @param gene_list List of character vectors of gene IDs.
#' @param rDEI_threshold rDEI threshold for filtering
#' @param upstream Number of base pairs upstream of the transcription start site
#' @param downstream Number of base pairs downstream of the transcription start site
#' @return Nothing returned
#' @export

genes_to_promoter_Biostring <- function(
  gene_list,
  rDEI_threshold = 1,
  upstream = 1000,
  downstream = 500
){
  ## Filter genes by rDEI if necessary
  gene_list <- lapply(gene_list, function(genes){filter_genes_by_rDEI(geneID = genes, rDEI_threshold = rDEI_threshold)})

  ## Convert genes to transcripts
  transcript_list <- lapply(gene_list, function(x){paste0(x, ".1")})

  ## Get the actual sequences
  Biostring_list <- transcript_list %>%
    lapply(function(x){get_promoters(x, upstream = upstream, downstream = downstream)})

  return(Biostring_list)
}



