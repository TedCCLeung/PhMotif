#' Function to clean up .meme motif files. Remove duplicate names if there is any. Remove duplicate motifs if there is any.
#'
#' @importFrom BiocGenerics %in%
#' @importFrom BiocGenerics match
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @import BiocGenerics
#' @importFrom rlang .data
#'
#' @param qry_fa Character vector of file paths of .fa files to be mapped on.
#' @param bkg_fa Character vector of file paths of .fa files to be used as background.
#' @param qry_gff3 Character vector of file paths of .gff3 files to be mapped on.
#' @param bkg_gff3 Character vector of file paths of .gff3 files to be used as background.
#' @param alternative Character. Alternative hypothesis of the Fisher's test from the stats package.
#' @param motif_max_p Numerical. What value of p value to cap if very low p values were encountered.
#'
#' @return Data frame of the p value and fold enrichment.
#' @export

enrich_motifs_fisher <- function(
  qry_gff3,
  bkg_gff3,
  qry_fa,
  bkg_fa,
  alternative = "two.sided",
  motif_max_p = 1e-5
){

  ## Process the gff3 files ---

  Grange_qry <- rtracklayer::import(qry_gff3)
  Grange_bkg <- rtracklayer::import(bkg_gff3)
  Grange_bkg <- Grange_bkg[!Grange_bkg %in% Grange_qry]   ## Remove from background those in the query (for Fisher's)

  qry_df <- Grange_qry %>%
    as.data.frame() %>%
    dplyr::mutate(pvalue = as.numeric(.data$pvalue)) %>%
    dplyr::filter(as.numeric(.data$pvalue) < motif_max_p)

  bkg_df <- Grange_bkg %>%
    as.data.frame() %>%
    dplyr::mutate(pvalue = as.numeric(.data$pvalue)) %>%
    dplyr::filter(as.numeric(.data$pvalue) < motif_max_p)

  ## Remove the motif ID column and then check for duplicates
  #qry_duplicate_filter <- duplicated(qry_df[, !c(colnames(qry_df) %in% c("motif.i"))])
  #bkg_duplicate_filter <- duplicated(bkg_df[, !c(colnames(bkg_df) %in% c("motif.i"))])

  ## Process the fa files ---
  qry_biostr <- Biostrings::readDNAStringSet(qry_fa)
  bkg_biostr <- Biostrings::readDNAStringSet(bkg_fa)
  bkg_biostr <- bkg_biostr[!bkg_biostr %in% qry_biostr]

  ##############################################

  motif_names <- sort(unique(qry_df$motif))

  df_out <- do.call("rbind", lapply(motif_names, enrich_motifs_fisher_subfunc, qry_biostr, bkg_biostr, qry_df, bkg_df, alternative)) %>%
    as.data.frame()
  df_out <- cbind(motif_names, df_out)

  return(df_out)
}




enrich_motifs_fisher_subfunc <- function(
  motif,
  qry_biostr,
  bkg_biostr,
  qry_df,
  bkg_df,
  alternative = "two.sided",
  RC = TRUE,
  return_cont_table = FALSE
){

  ##############################################
  ## Select the data for the specific motif
  qry_df_sel <- qry_df[qry_df$motif == motif, ]
  bkg_df_sel <- bkg_df[bkg_df$motif == motif, ]

  motif_width <- mean(qry_df_sel$width)

  qry_max_possible_hits <- (mean(width(qry_biostr)) - motif_width + 1) * length(qry_biostr)
  bkg_max_possible_hits <- (mean(width(bkg_biostr)) - motif_width + 1) * length(bkg_biostr)

  if (RC) {
    qry_max_possible_hits <- qry_max_possible_hits * 2
    bkg_max_possible_hits <- bkg_max_possible_hits * 2
  }

  bkg_norm_factor <- qry_max_possible_hits/bkg_max_possible_hits

  qry_hit <- nrow(qry_df_sel)
  qry_no_hit <- qry_max_possible_hits - qry_hit
  bkg_hit <- nrow(bkg_df_sel)
  bkg_no_hit <- (bkg_max_possible_hits - bkg_hit)

  ##############################################
  ## Build the contingency table of the Fisher's exact test
  contingency_mat <- matrix(c(
    qry_hit, qry_no_hit,
    bkg_hit * bkg_norm_factor, bkg_no_hit * bkg_norm_factor
  ),
  nrow = 2,
  byrow = TRUE
  )
  ##############################################
  ## Fisher's exact test
  fisher_out <- stats::fisher.test(round(contingency_mat), alternative = alternative)
  fold_enrichment <- (qry_hit/qry_no_hit)/((bkg_hit+qry_hit)/(bkg_no_hit+qry_no_hit))

  if (return_cont_table) {
    contingency_mat_ <- matrix(c(qry_hit, qry_no_hit, bkg_hit, bkg_no_hit), ncol = 2, byrow = TRUE)
    return(contingency_mat_)
  } else {
    return(c("pvalue" = fisher_out$p.value, "fold_enrichment" = fold_enrichment))
  }
}


