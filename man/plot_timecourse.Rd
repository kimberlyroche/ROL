% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_timecourse}
\alias{plot_timecourse}
\title{Render bar plots for an individual host's series}
\usage{
plot_timecourse(
  data,
  host = NULL,
  gapped = FALSE,
  legend = TRUE,
  selected_samples = NULL,
  show_plot = FALSE
)
}
\arguments{
\item{data}{a phyloseq object}

\item{host}{host short name (e.g. ACA)}

\item{gapped}{render true sample spacing (approximately)}

\item{legend}{render legend}

\item{selected_samples}{list of sample IDs to highlight (if visualizing a selection of a host's samples)}

\item{show_plot}{show() plot in addition to rendering it to a file}
}
\description{
Render bar plots for an individual host's series
}
\examples{
plot_timecourse(data, host="ACA", gapped=FALSE, legend=TRUE, selected_samples=NULL)
}
