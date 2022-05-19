#' Find similar motifs against the CIS-BP motif database
#'
#' @param query_motifs Universalmotif objects
#'
#' @return Dataframe object
#'
#' @export

motif_to_similar_CISBP <- function(
  query_motifs
){

  ## Get query motif names
  query_motif_names <- query_motifs %>% universalmotif::to_df() %>% magrittr::extract2("name")

  ## Motif comparison
  motif_sim_mat <- universalmotif::compare_motifs(
    c(CISBP_motifset, query_motifs)
  )

  ## Remove the comparisons of query motifs against other query motifs
  filtered_sim_mat <- motif_sim_mat[!(rownames(motif_sim_mat) %in% query_motif_names), query_motif_names]

  ## Summarize the result into a dataframe
  df_sim_result <- data.frame(

    query_motif = query_motif_names,

    sim_motif = filtered_sim_mat %>%
      apply(2, function(x){names(sort(x[!is.na(x)], decreasing = TRUE)[1])}),

    sim = filtered_sim_mat %>%
      apply(2, function(x){sort(x[!is.na(x)], decreasing = TRUE)[1]})
  )

  return(df_sim_result)
}
