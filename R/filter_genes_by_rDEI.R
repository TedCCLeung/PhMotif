#' Function to filter gene IDs (TAIR IDs) based on daily expression integral
#'
#' @param geneID Character vector of genes in TAIR ID
#' @param rDEI_threshold rDEI threshold
#' @param rDEI_file File of rDEI values
#' @return Character vector of genes that pass the rDEI threshold
#' @export

filter_genes_by_rDEI <- function(
  geneID,
  rDEI_threshold = 1,
  rDEI_file = NULL
){

  if (is.null(rDEI_file)){df_rDEI <- rDEI_summary
  } else {df_rDEI <- utils::read.csv(rDEI_file)}

  rDEI_cols <- colnames(df_rDEI)[startsWith(colnames(df_rDEI), "log2_rDEI_")]
  max_rDEI <- suppressWarnings(apply(abs(df_rDEI[, rDEI_cols]), 1, function(m){max(m, na.rm = TRUE)}))
  max_rDEI[is.infinite(max_rDEI) | is.na(max_rDEI)] <- 0
  high_rDEI_genes <- df_rDEI[max_rDEI > log2(rDEI_threshold), "geneID"]
  genes <- geneID[geneID %in% high_rDEI_genes]
  return(genes)
}
