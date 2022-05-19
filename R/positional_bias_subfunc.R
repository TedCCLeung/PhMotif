#' Analyse positional bias of motifs in a set of DNA sequences
#'
#' @param motifs Universalmotif object or name of a .meme file with motifs
#' @param list_of_seqs List of Biostrings objects or list of file names of .fa file for DNA sequences
#' @param seq_names Names of the individual items in the list_of_seqs
#' @param motif_name Character. Name of the motif
#' @param p_value Numerical. p-value for calling a motif peak
#' @param reverse_complement Logical. Whether the reverse complement is used to look for motifs
#' @param bandwidth Bandwidth of the universalmotif package
#'
#' @return Data frame
#'

positional_bias_subfunc <- function(
  motifs,
  list_of_seqs,
  seq_names = NULL,
  motif_name = "NA",
  p_value = 1e-06,
  reverse_complement = TRUE,
  bandwidth = NULL
){

  result_list <- lapply(1:length(list_of_seqs), function(k){

    ## Get the sequence and the names
    seqs <- list_of_seqs[[k]]
    seqs_name <- ifelse(is.null(seq_names), as.character(k), seq_names[k])

    hits <- universalmotif::scan_sequences(motifs, seqs, RC = reverse_complement)

    if (is.null(bandwidth)){
      res <- universalmotif::motif_peaks(
        hits$start,
        max.p = p_value
      )
    } else {
      res <- universalmotif::motif_peaks(
        hits$start,
        max.p = p_value,
        bandwidth = bandwidth
      )
    }

    ## the motif_peaks function returns 2 items if a peak is found, and just the plot if no significant peak is found
    if (length(res) == 2){output_table <- res$Plot$data} else {output_table <- res$data}
    #output_table <- ifelse(length(res) == 2, res$Plot$data, res$data)
    output_table$significant <- rep(length(res) == 2, nrow(output_table))
    output_table$sequence <- rep(seqs_name, nrow(output_table))
    output_table$motif <- rep(motif_name, nrow(output_table))

    return(output_table)
  })

  ## combine the distribution dfs into a single one for plotting
  df <- do.call("rbind", result_list)

  return(df)
}
