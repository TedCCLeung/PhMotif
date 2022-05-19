#' Perform hierarchical clustering of motifs
#'
#' @param motif_family Character.
#' @param custom_TF_family_scheme a custom naming scheme for the TF families as illustrated in df_TFfamily_colors
#'
#' @return List of ComplexHeatmap objects
#'
#' @export

heatmap_TF_family <- function(
  motif_family,
  custom_TF_family_scheme = NULL
){

  ## COLOR SCHEME --------------------------------

  ## Check whether there is a custom color and TF family labeling scheme is supplied
  if (is.null(custom_TF_family_scheme)){
    color_scheme <- df_TFfamily_colors$color
    names(color_scheme) <- df_TFfamily_colors$label
    motif_group <- df_TFfamily_colors[match(motif_family, df_TFfamily_colors$family), "label"]
  } else {
    custom_scheme <- custom_TF_family_scheme
    color_scheme <- custom_scheme$color
    names(color_scheme) <- custom_scheme$label
    motif_group <- custom_scheme[match(motif_family, custom_scheme$family), "label"]
  }

  ## HEATMAP BAR --------------------------------

  bar <- ComplexHeatmap::Heatmap(
    ## Input
    motif_group,
    ## Color
    col = color_scheme,
    name = "Family"
  )

  ## RETURN RESULTS --------------------------------
  return(bar)
}

