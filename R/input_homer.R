input_homer <- function(
  dir
){
  motif_files <- list.files(path = dir, recursive = TRUE, pattern = "knownResults.txt", all.files = TRUE, full.names = TRUE)
  cond_names <- motif_files %>% strsplit(split = "/") %>% sapply(function(x){magrittr::extract2(x, length(x)-1)})

  df_result <- Reduce(rbind, Map(read_homer_output, motif_files, cond_names))

  return(df_result)
}
