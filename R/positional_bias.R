#' Analyse positional bias of motifs in a set of DNA sequences
#'
#' @importFrom magrittr %>%
#'
#' @param motifs Universalmotif object or name of a .meme file with motifs
#' @param seqs List of Biostrings objects or list of file names of .fa file for DNA sequences
#' @param separate_motifs Logical. Whether individual motif should be analyzed separately
#' @param seq_names Names of the individual items in the list_of_seqs
#' @param pvalue Numerical. p-value for calling a motif peak
#' @param reverse_complement Logical. Whether the reverse complement is used to look for motifs
#' @param bandwidth Bandwidth of the universalmotif package
#'
#' @return Data frame
#'
#' @export

positional_bias <- function(
  motifs,
  seqs,
  separate_motifs = TRUE,
  seq_names = NULL,
  pvalue = 1e-06,
  reverse_complement = TRUE,
  bandwidth = NULL
){

  ## Get the sequences (in case they are input as list of file names)
  if (is.character(seqs)){
    list_of_seqs <- lapply(seqs, function(x){y <- Biostrings::readDNAStringSet(x); y[lengths(y) == get_mode(lengths(y))]})
  } else {
    list_of_seqs <- seqs
  }

  ## In case the motifs is input as a file name
  if (is.character(motifs)){
    motifs <- universalmotif::read_meme(motifs)
  } else {
    motifs <- motifs
  }

  ## If the seq_names argument is not filled, and the list_of_seqs object has fitting names, use them as names
  if (is.null(seq_names)){
    if (!is.null(names(list_of_seqs))){
      if (length(names(list_of_seqs)) == length(list_of_seqs)){
        seq_names <- names(list_of_seqs)
      }
    }
  }

  ## In the case of treating motifs separately -------------------------------
  if (separate_motifs){

    motif_result <- lapply(motifs, function(motif){

      position_result <- positional_bias_subfunc(
        motifs = motif,
        list_of_seqs = list_of_seqs,
        p_value = pvalue,
        seq_names = seq_names,
        motif_name = motif@name,
        reverse_complement = reverse_complement,
        bandwidth = bandwidth
      )
      return(position_result)
    })

    position_result <- do.call("rbind", motif_result)
    return(position_result)

    ## In the case of treating motifs together -------------------------------
  } else {

    position_result <- positional_bias_subfunc(
      motifs = motifs,
      list_of_seqs = list_of_seqs,
      p_value = pvalue,
      seq_names = seq_names,
      motif_name = "NA",
      reverse_complement = reverse_complement,
      bandwidth = bandwidth
    )
    return(position_result)
  }
}


