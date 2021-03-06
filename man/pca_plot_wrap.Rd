% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_wrappers.R
\name{pca_plot_wrap}
\alias{pca_plot_wrap}
\title{PCA plot wrapper}
\usage{
pca_plot_wrap(
  count_set = NULL,
  norm_method = "housekeeping_scaled",
  title = NA,
  comp1 = 1,
  comp2 = 2,
  legend = TRUE,
  label = TRUE,
  colors = NA
)
}
\arguments{
\item{count_set}{A count_set of mRNA or miRNA counts generated using count_set}

\item{norm_method}{Name of the count_set assay to be plotted. Options are
"counts", "background_corrected","positive_control_scaled",
"housekeeping_scaled", "geNorm_housekeeping", "all_endogenous_scaled", "loess", "vsn", "quantile", "ruvIII".
Default = "housekeeping_scaled"}

\item{title}{Plot title}

\item{comp1}{PCA compenent to plot on the x axis. Default = Comp 1}

\item{comp2}{PCA compenent to plot on the y axis Default = Comp 1}

\item{legend}{Option for plot legend. Default = TRUE}

\item{label}{Option to label points with sample id. Default= TRUE}

\item{colors}{Option to define colour scheme. Default = NA}
}
\value{
pca plots
}
\description{
PCA plot wrapper
}
