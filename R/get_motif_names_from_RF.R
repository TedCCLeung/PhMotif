get_motif_names_from_RF <- function(
  file = "/Users/TedCCLeung/Documents/Projects/Photoperiod/motif_importance.txt",
  cutoff = 0.001
){

  df <- utils::read.table(file, header = TRUE, sep = "\t")
  df <- df[order(df$importance, decreasing = TRUE), ]
  motif_names <- df[df$importance > cutoff, "motif"]

  return(motif_names)
}
