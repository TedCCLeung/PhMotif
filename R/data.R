#' Motif set formed from de novo motif analysis with WEEDER, HOMER and STREME.
#'
#' @format MEME format
"denovo_motifset1"

#' Motif set from the JASPAR Plants database.
#'
#' @format MEME format
"JASPAR_motifset"

#' TAIR10 gene models from Arabidopsis_thaliana.TAIR10.47.gtf
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{EQ_DEI}{daily expression integral at equinox photoperiod}
#'   \item{LD_DEI}{daily expression integral at long day photoperiod}
#'   \item{SD_DEI}{daily expression integral at short day photoperiod}
#'   \item{geneID}{}
#'   \item{log2_rDEI_LDEQ}{}
#'   \item{log2_rDEI_SDEQ}{}
#'   \item{log2_rDEI_SDLD}{}
#'   \item{rDEI_LDEQ}{}
#'   \item{rDEI_SDEQ}{}
#'   \item{rDEI_SDLD}{}
#'
#' }
"rDEI_summary"

#' Grange object from Arabidopsis_thaliana.TAIR10.47.gtf
#'
#' @source \url{https://plants.ensembl.org/Arabidopsis_thaliana/Info/Index}
"TAIR10_Grange"


#' Motif set from the CIS-BP database.
#'
#' @format MEME format
#' @source \url{http://cisbp.ccbr.utoronto.ca/index.php}
"CISBP_motifset"

#' Dataframe of motif families from the CIS-BP database.
#'
#' @format Data frame
#' #' \describe{
#'   \item{family}{TF family as determined by CISBP}
#'   \item{Freq}{Frequency of transcription factor motifs in CISBP}
#'   \item{label}{Name of TF family that is to be displayed}
#'   \item{color}{RGB code of the family label}
#' }
"df_TFfamily_colors"

#' Dataframe of families of motifs from the CIS-BP database.
#'
#' @format Data frame
#' #' \describe{
#'   \item{motif}{Name of the motifs from CIS-BP database}
#'   \item{family}{TF family as determined by CISBP}
#' }
"CISBP_families"

#' Photoperiodic clusters
#' @format List of list of characters
"photoperiodic_clusters"
