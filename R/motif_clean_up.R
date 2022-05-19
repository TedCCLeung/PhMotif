#' Function to clean up .meme motif files. Remove duplicate names if there is any. Remove duplicate motifs if there is any.
#'
#' @param input_files Character vector of file paths of .meme files to be cleaned up.
#' @param keep_duplicates Whether to keep motifs with the same positional weight matrices.
#' @param keep_same_name_motifs Whether to keep motifs with the same name.
#' @param output_file Whether to write a .meme motif file.
#' @return A universalmotif object or nothing if output_file=TRUE.
#' @export

motif_clean_up <- function(
  input_files,
  keep_duplicates = FALSE,
  keep_same_name_motifs = TRUE,
  output_file = NULL
){

  ## Convert the motifs into a data frame for easy manipulation
  raw_motifs <- sapply(input_files, universalmotif::read_meme) %>% unlist()

  if (!keep_duplicates){
    df_motifs <- universalmotif::merge_similar(raw_motifs, threshold = 1.00) %>% universalmotif::to_df()
    print("Number of motifs removed:")
    print(length(raw_motifs) - nrow(df_motifs))
    ## Simplify the altnames column
    df_motifs$altname <- df_motifs$altname %>% sapply(function(x){return(utils::head(strsplit(x, split = "/")[[1]], 1))})
  } else {
    df_motifs <- raw_motifs %>% universalmotif::to_df()
  }

  ## If there duplicated motif names
  if (sum(duplicated(df_motifs$name)) > 0 | keep_same_name_motifs == FALSE){
    if (sum(duplicated(df_motifs$altname)) == 0) {
      print("There are no duplicated altnames. Altnames will be switched with names.")
      df_motifs %<>% dplyr::rename("altname" = "name", "name" = "altname")
    } else {
      print("There are duplicated altnames. Names will be modified")
      df_motifs$name <- sub("\\.", "_", make.names(df_motifs$name))
    }
  }

  motifs_out <- universalmotif::to_list(df_motifs, extrainfo = FALSE)

  if (is.null(output_file)){
    return(motifs_out)
  } else {
    universalmotif::write_meme(motifs_out, output_file)
  }

}
