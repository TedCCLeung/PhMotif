#' Function to return Biostring object with input of gene models
#'
#' @param upstream Numerical. Number of base pair upstream of the transcription start site.
#' @param downstream Numerical. Number of base pair downstream of the transcription start site.
#' @param gene_models Character. Name of TAIR10 transcripts.
#'
#' @export

get_promoters <- function(
  gene_models,
  upstream,
  downstream
){

  ## 01. Get the location of promoter
  ## Ideally I do not want to make a TxDb object everytime I get promoters
  ## But somehow the sqlite database doesn't work when but inside a package due to null pointer issues
  txdb <- make_TAIR10_TxDb()
  valid.genes <- GenomicFeatures::transcripts(txdb)$tx_name
  gene_models <- gene_models[gene_models %in% valid.genes]
  print(paste0("invalid genes: ", as.character(length(gene_models[!(gene_models %in% valid.genes)]))))
  pr.ranges <- GenomicFeatures::promoters(txdb, upstream = upstream, downstream = downstream)[gene_models,]
  ## Add the length of chromosomes in bp manually
  GenomeInfoDb::seqlengths(pr.ranges) <- c(30427671, 19698289, 23459830, 18585056, 26975502, NA, NA)
  pr.ranges <- GenomicRanges::trim(pr.ranges)
  GenomeInfoDb::seqlevelsStyle(pr.ranges) <- "TAIR9"

  ## 02. Get the genome sequence (TAIR9 and TAIR10 are identical)
  genome <- BSgenome.Athaliana.TAIR.TAIR9::BSgenome.Athaliana.TAIR.TAIR9

  ## 03. Get the sequence and output to a fasta file
  return(Biostrings::getSeq(genome, pr.ranges))
}

get_fiveUTR <- function(
  gene_models
){
  ## 01. Get the location of promoter
  txdb <- make_TAIR10_TxDb()
  valid.genes <- GenomicFeatures::transcripts(txdb)$tx_name
  gene_models <- gene_models[gene_models %in% valid.genes]
  print(paste0("invalid genes: ", as.character(length(gene_models[!(gene_models %in% valid.genes)]))))
  five.utr <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
  genes.with.five <- names(five.utr)
  gene_models <- gene_models[gene_models %in% genes.with.five]
  print(paste0("genes with unknown 5'UTR: ", as.character(gene_models[!(gene_models %in% genes.with.five)])))
  five.utr <- five.utr[gene_models]
  GenomeInfoDb::seqlevelsStyle(five.utr) <- "TAIR9"
  GenomeInfoDb::seqlengths(five.utr) <- c(30427671, 19698289, 23459830, 18585056, 26975502, NA, NA)
  five.utr <- GenomicRanges::trim(five.utr)
  five.utr <- unlist(five.utr)

  ## 02. Get the genome sequence (TAIR9 and TAIR10 are identical)
  genome <- BSgenome.Athaliana.TAIR.TAIR9::BSgenome.Athaliana.TAIR.TAIR9

  ## 03. Get the sequence and output to a fasta file
  return(Biostrings::getSeq(genome, five.utr))
}

get_prom_five <- function(
  gene_models,
  prom.upstream
){

  prom <- get_promoters(
    upstream = prom.upstream,
    downstream = 0,
    gene_models = gene_models
  )

  five <- get_fiveUTR(
    gene_models = gene_models
  )

  res <- Biostrings::DNAStringSet(x = "A")
  for (k in 1:length(gene_models)){
    gene <- gene_models[k]
    if ((gene %in% names(prom)) & (gene %in% names(five))){
      gene.res <- Biostrings::xscat(prom[gene], five[gene])
      names(gene.res) <- gene
      res <- c(res, gene.res)
    } else if ((gene %in% names(prom))){
      res <- c(res, prom[gene])
    }
  }

  return(res[2:length(res)])
}
