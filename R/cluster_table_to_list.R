#' Convert a cluster table to a list of list object
#'
#' @importFrom magrittr %>%
#'
#' @param cluster_table_file Character. Full path to output from the package PhotoperiodCluster
#'
#' @return List of list


cluster_table_to_list <- function(
  cluster_table_file
){

  df_input <- utils::read.table(cluster_table_file, header = TRUE, sep = "\t")

  big_clusters <- df_input$static %>% unique()

  get_sub_clusters <- function(y){df_input[df_input$sub == y, "geneID"]}

  final_res <- lapply(big_clusters, function(x){
    big_df <- df_input[df_input$static == x, ]
    sub_res <- lapply(unique(big_df$sub), get_sub_clusters)
    names(sub_res) <- unique(big_df$sub)
    return(sub_res)
  })

  names(final_res) <- big_clusters %>% addLeadingZeros()

  return(final_res)
}
