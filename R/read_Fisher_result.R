#' Retrieve results from Fisher's Exact Test
#'
#' @param in_dir Directory with the outputs of processFisherData
#' @param max_pval_cap Numerical. The maximum p value to display on a heatmap
#'
#' @return A list of mat_logp, mat_logfdr, mat_FC and motif_names
#'
#' @export

read_Fisher_result <- function(
  in_dir,
  max_pval_cap = 100
){

  ## PROCESS THE P VALUE MATRIX -----------------------------
  pval_file <- paste0(in_dir, "/pvalue.tsv")
  mat_pval <- pval_file %>%
    utils::read.table(sep = "\t", header = TRUE) %>%
    tibble::column_to_rownames(var = "motif_names") %>%
    as.matrix()
  mat_logp <- log10(mat_pval) * -1
  mat_logp[is.infinite(mat_logp)] <- max_pval_cap

  ## PROCESS THE FDR MATRIX -----------------------------
  fdr_file <- paste0(in_dir, "/fdr.tsv")
  mat_fdr <- fdr_file %>%
    utils::read.table(sep = "\t", header = TRUE) %>%
    tibble::column_to_rownames(var = "motif_names") %>%
    as.matrix()
  mat_logfdr <- log10(mat_fdr) * -1
  mat_logfdr[is.infinite(mat_logfdr)] <- max_pval_cap

  ## PROCESS THE FC -----------------------------
  fc_file <- paste0(in_dir, "/fold_enrichment.tsv")
  mat_FC <- fc_file %>%
    utils::read.table(sep = "\t", header = TRUE) %>%
    tibble::column_to_rownames(var = "motif_names") %>%
    as.matrix() %>%
    log2()

  ## GET UNIVERSALMOTIF OBJECTS -----------------------------
  motif_names <- rownames(mat_pval)

  result <- list(
    "mat_logPVAL" = mat_logp,
    "mat_logFDR" = mat_logfdr,
    "mat_FC" = mat_FC,
    "motif_names" = motif_names
  )

  return(result)
}
