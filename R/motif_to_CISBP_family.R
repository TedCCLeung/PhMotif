#' Search which family the most motif belongs to in CISBP
#'
#' @param motifs Universalmotif objects
#' @param similarity_cutoff Pearson's similarity cutoff for determining if a motif has a counterpart in CISBP
#' @param only_CISBP_motifs Whether the input consists of only motifs from CISBP and comparison is unnecessary
#'
#' @return Character
#'
#' @export

motif_to_CISBP_family <- function(
  motifs,
  similarity_cutoff = 0.5,
  only_CISBP_motifs = FALSE
){

  if (only_CISBP_motifs){

    ## This skips motif comparison if the input motifs are all from CISBP

    ## Get the name of input motifs
    query_motif_names <- motifs %>% universalmotif::to_df() %>% magrittr::extract2("name")
    ## Search for the family using the CISBP_families dataframe
    result_families <- CISBP_families$family[match(query_motif_names, CISBP_families$motif)]
    return(result_families)

  } else {

    ## Search against the CISBP database
    search_result <- motif_to_similar_CISBP(query_motifs = motifs)

    ## Get the name of input motifs
    query_motif_names <- motifs %>% universalmotif::to_df() %>% magrittr::extract2("name")

    ## Get the family of the most similar CISBP motif
    most_similar_motifs <- search_result[query_motif_names, "sim_motif"]
    ## Remove the potential for getting duplication of CISBP in the similarity matrix
    most_similar_motifs <- gsub(" \\[\\s*(.*?)\\s*\\]", "", most_similar_motifs)

    ## Search for the family using the CISBP_families dataframe
    result_families <- CISBP_families$family[match(most_similar_motifs, CISBP_families$motif)]

    ## Convert the resulted family to "Unknown" if the similarity value is below the cutoff
    sim_filter <- search_result[query_motif_names, "sim"] < similarity_cutoff
    result_families[sim_filter] <- "NoMatch"

    return(result_families)

  }

}
