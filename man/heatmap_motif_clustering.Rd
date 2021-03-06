% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap_motif_clustering.R
\name{heatmap_motif_clustering}
\alias{heatmap_motif_clustering}
\title{Plot the result of hierarchical clustering of motifs}
\usage{
heatmap_motif_clustering(
  input_matrix,
  hclust_out,
  cutree_out,
  motifset,
  colors = c("purple", "white", "orange"),
  color_range = c(-2, 0, 2),
  scale_name = " ",
  panel_size = grid::unit(0.75, "cm"),
  panel_gap = grid::unit(0, "cm"),
  panel_width = grid::unit(4, "cm"),
  panel_link_width = grid::unit(0.6, "cm"),
  remove_zoom = FALSE
)
}
\arguments{
\item{input_matrix}{Numerical matrix with motif names as row names and sample names as column names.}

\item{hclust_out}{result from hclust() of the motifs}

\item{cutree_out}{result from cutree()}

\item{motifset}{Universalmotif object}

\item{colors}{Character vector of colors for the gradient.}

\item{color_range}{Numerical vector that represents the values the number correspond to. Must be of the same length as colors.}

\item{scale_name}{Title of the gradient scale legend.}

\item{panel_size}{Size of the panel in ggplot::grid::unit().}

\item{panel_gap}{Gap of the panel in ggplot::grid::unit().}

\item{panel_width}{Width of the panel in ggplot::grid::unit().}

\item{panel_link_width}{Link width of the panel in ggplot::grid::unit().}

\item{remove_zoom}{whether to remove zoom box}
}
\value{
ComplexHeatmap object
}
\description{
Plot the result of hierarchical clustering of motifs
}
