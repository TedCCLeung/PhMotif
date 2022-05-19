## -------------------------------------------
#' Processing of results from fisherEnrichMotif in a directory
#' @param dir Character. Input directory. Often the ouput directory of fisherEnrichMotif_inDir.
#' @param max_pval Numerical. p-value threshold above which the data will not be incoporated in the processed tables.
#' @param outdir Character. Output directory.
#' @return None
#' @export
enrich_result_summary <- function(
  dir,
  max_pval = 0.05,
  outdir
){
  ## Retrieve the enrichment data
  FisherData_out <- pivotFisherData(input_dir = dir, max_pval = max_pval)
  df_pval <- FisherData_out$pval
  df_fdr <- FisherData_out$BH
  df_FC <- FisherData_out$FC
  utils::write.table(df_pval, file = paste0(outdir, "pvalue.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(df_fdr, file = paste0(outdir, "fdr.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(df_FC, file = paste0(outdir, "fold_enrichment.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}



## -------------------------------------------
## Function to pivot Fisher data in a directory
pivotFisherData <- function(
  input_dir,
  max_pval = NULL
){

  df_FC <- pivotDFinDir(
    input_dir = input_dir,
    ID_col_name = "motif_names",
    value_col_name = "fold_enrichment"
  )
  df_FC[is.na(df_FC)] <- 1

  df_pval <- pivotDFinDir(
    input_dir = input_dir,
    ID_col_name = "motif_names",
    value_col_name = "pvalue"
  )
  df_pval[is.na(df_pval)] <- 1

  df_BH_ <- df_pval[, 2:ncol(df_pval)] %>% as.matrix() %>% as.numeric() %>% stats::p.adjust(method = "BH") %>%
    matrix(ncol = ncol(df_pval)-1, byrow = FALSE) %>% as.data.frame()
  df_BH <- cbind(df_pval$motif_names, df_BH_)
  colnames(df_BH) <- colnames(df_pval)

  if (!is.null(max_pval)){
    pval_filter <- apply(df_pval[, 2:ncol(df_pval)], 1, min) < max_pval
    BH_filter <- rowSums(df_BH[, 2:ncol(df_BH)]) != ncol(df_BH)-1
    df_FC <- df_FC[pval_filter & BH_filter, ]
    df_pval <- df_pval[pval_filter & BH_filter, ]
    df_BH <- df_BH[pval_filter & BH_filter, ]
  }

  return(list("FC" = df_FC, "pval" = df_pval, "BH" = df_BH))
}





