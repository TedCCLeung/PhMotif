#' Function to generate a TxDb object from the GenomicFeatures object for getting gene promoters.
#'
#' @return Returns a TxDb object (from the GenomicFeatures package)
#' @export

make_TAIR10_TxDb <- function(
){
  tair10 <- GenomicFeatures::makeTxDbFromGRanges(TAIR10_Grange)
  return(tair10)
}
