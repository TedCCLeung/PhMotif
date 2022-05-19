## -------------------------------------------
## Function used in test_positional_bias.R
get_mode <- function(x) {
  y <- unique(x)
  y[which.max(tabulate(match(x, y)))]
}

## -------------------------------------------
#' Read a .txt file line by line and return a character vector.
#' @param x Character. Name of a .txt file
#' @return Character.
#' @export
read_genes <- function(x){y <- readLines(file(x)); close(file(x)); return(sort(y))}

## -------------------------------------------
#' Add leading zeros to a numerical vector and convert it to string.
#' @param num_vec Numerical. Numerical vector.
#' @return Character.
#' @export
addLeadingZeros <- function(
  num_vec
){
  digit_number <- floor(max(log10(num_vec)+1))
  return(stringr::str_pad(as.character(num_vec), digit_number, pad = "0"))
}

## -------------------------------------------
#' Concatenate a list of FASTA files
#' @importFrom magrittr %>%
#' @param input_FASTA Character. Input FASTA file(s).
#' @param output_FASTA Character. Output FASTA file.
#' @return None
#' @export
concatenate_FASTA <- function(
  input_FASTA,
  output_FASTA
){
  do.call("c", lapply(input_FASTA, Biostrings::readDNAStringSet)) %>%
    Biostrings::writeXStringSet(filepath = output_FASTA)
}


## -------------------------------------------
#' Concatenate a list of GFF3 files
#' @importFrom magrittr %>%
#' @param input_gff3 Character. Input GFF3 file(s).
#' @param output_gff3 Character. Output GFF3 file.
#' @return None
#' @export
concatenate_gff3 <- function(
  input_gff3,
  output_gff3
){
  suppressWarnings(do.call("c", lapply(input_gff3, rtracklayer::import)) %>% rtracklayer::export(output_gff3, format = "gff3"))
}

## -------------------------------------------
#' Merge all the tabular text files in a directory and perform an optional pivot
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @param input_dir Character. Name of directory.
#' @param ID_col_name Character. Name of the column that contains the ID to be used as pivot.
#' @param value_col_name Character. Name of the column that contains the value to be used as pivot.
#' @param sep Character. Separator of the tabular data.
#' @param skip_pivot Logical. Whether the pivot should be performed.
#' @return Data frame
#' @export

pivotDFinDir <- function(
  input_dir,
  ID_col_name,
  value_col_name,
  sep = "\t",
  skip_pivot = FALSE
){
  group_names <- list.files(input_dir) %>% sapply(function(x){utils::head(strsplit(x, split = "\\.")[[1]], 1)})
  df_list <-  list.files(input_dir, full.names = TRUE) %>%
    lapply(function(x){utils::read.table(x, sep = sep, header = TRUE)[, c(ID_col_name, value_col_name)]})
  df_all <- do.call("rbind", df_list)
  df_all$group <- rep(group_names, times = sapply(df_list, nrow))

  if (skip_pivot){
    return(df_all)
  } else {
    df_pivot <- df_all %>%
      dplyr::select(.data$group, dplyr::all_of(ID_col_name), dplyr::all_of(value_col_name)) %>%
      tidyr::pivot_wider(names_from = .data$group, values_from = dplyr::all_of(value_col_name))
    return(df_pivot)
  }
}


## --------

get_file_names <- function(
  gene_files
){
  file_names <- lapply(gene_files, function(x){
    f_name <- stringr::str_trim(utils::tail(strsplit(x, split = "/")[[1]], n = 1))
    out_name <- stringr::str_trim(paste0(utils::head(strsplit(f_name, split = "\\.")[[1]], n = 1), ".fa"))
    return(out_name)
  })
  return(file_names)
}


## -------
## Function to convert MEME file paths into universalmotifs
input_motifs <- function(
  motifs,
  out_file = NULL
){
  ## Check if the input is a universalmotif object of a path to a MEME file
  if (class(motifs[[1]]) == "universalmotif"){motifs_ <- motifs
  } else if (class(motifs[[1]]) == "character"){motifs_ <- universalmotif::read_meme(motifs)
  } else {stop("Wrong motif input!")
  }

  if (!is.null(out_file)){
    universalmotif::write_meme(motifs = motifs_, file = out_file, overwrite = TRUE)
  }

  return(motifs_)
}


read_homer_output <- function(
  x,
  sample
){
  df <- utils::read.delim(x)

  percent_target <- gsub('.{1}$', '', df[, 7]) %>% as.numeric()
  percent_background <- gsub('.{1}$', '', df[, 9]) %>% as.numeric()
  fold_change <- log2(percent_target/percent_background)

  pvalue <- exp(df[, 4] %>% as.numeric())
  fdr <- df[, 5] %>% as.numeric()
  name <- df[, 1]

  df_out <- data.frame(
    name,
    sample = rep(sample, length(name)),
    fold_change,
    pvalue,
    fdr
  )
  return(df_out)
}


input_homer <- function(
  dir
){
  motif_files <- list.files(path = dir, recursive = TRUE, pattern = "knownResults.txt", all.files = TRUE, full.names = TRUE)

  cond_names <- motif_files %>%
    strsplit(split = "/") %>%
    sapply(function(x){magrittr::extract2(x, length(x)-1)})

  df_result <- Reduce(rbind, Map(read_homer_output, motif_files, cond_names))

  return(df_result)
}

