% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cytoMerge.R
\name{freqplot}
\alias{freqplot}
\title{Heatmap for 2D histograms}
\usage{
freqplot(
  x,
  y,
  breaks = 200,
  na.rm = TRUE,
  palette = rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF",
    "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")),
  add_white = TRUE,
  ...
)
}
\arguments{
\item{x}{numerical vector}

\item{y}{numerical vector}

\item{breaks}{Number of bins}

\item{na.rm}{Whether to exclude NAs from x and y}

\item{palette}{color palette for the heatmap}

\item{add_white}{Boolean. Whether to add white at the beginning of the palette (for counts of 0 in the 2D histogram)}

\item{...}{passed to image()}
}
\description{
Heatmap for 2D histograms
}
