#' Perform hierarchical clustering of motifs
#'
#' @param motifs Universalmotif object.
#' @param compare_motifs_method The similarity metric to be used between motifs according to the Universalmotifs package
#' @param hclust_method Method to perform hierarchical clustering.
#' @param min.position.ic The minimum IC between positions according to the Universalmotifs package
#' @param min.mean.ic The minimum mean IC according to the Universalmotifs package
#' @param min.overlap The minimum overlap between motifs according to the Universalmotifs package
#' @param normalize.scores Whether to normalize scores according to the Universalmotifs package
#'
#' @return hclust result
#'
#' @export

hclust_motifs <- function(
  motifs,
  compare_motifs_method = "PCC",
  hclust_method = "complete",
  min.position.ic = 0,
  min.mean.ic = 0.25,
  min.overlap = 6,
  normalize.scores = TRUE
){
  mat_motifs <- 1 - universalmotif::compare_motifs(
    motifs,
    method = compare_motifs_method,
    normalise.scores = normalize.scores,
    min.position.ic = min.position.ic,
    min.mean.ic = min.mean.ic,
    min.overlap = min.overlap
  )
  mat_motifs[is.na(mat_motifs)] <- max(mat_motifs, na.rm = TRUE)
  hclust_out <- stats::hclust(stats::as.dist(mat_motifs), method = hclust_method)
  return(hclust_out)
}
