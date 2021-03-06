% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_fitting.R
\name{perform_cv}
\alias{perform_cv}
\title{Perform k-fold cross-validation for a given host and kernel configuration choice}
\usage{
perform_cv(
  data,
  host,
  taxa_covariance,
  sample_covariance,
  tax_level = "ASV",
  holdout_proportion = 0.2
)
}
\arguments{
\item{data}{a phyloseq object}

\item{host}{host short name (e.g. ACA)}

\item{sample_covariance}{composite kernel function}

\item{tax_level}{taxonomic level at which to agglomerate data}

\item{holdout_proportion}{proportion of host's sample to use as a test set}

\item{rho}{scalar bandwidth parameter}

\item{iterations}{number of CV iterations to perform; this should be <= n_samples for this host}
}
\description{
Perform k-fold cross-validation for a given host and kernel configuration choice
}
\details{
average error (currently RMSE of log counts)
}
\examples{
tax_level <- "ASV"
data <- load_data(tax_level = tax_level)
sample_covariance <- get_Gamma(kernel_scale = 2, proportions = c(1, 0, 0), min_correlation = 0.1, days_to_baseline = 90)
taxa_covariance <- get_Xi(phyloseq::ntaxa(data), total_variance = 1)
perform_cv(data, host = "GAB", taxa_covariance = taxa_covariance, sample_covariance = sample_covariance, tax_level = tax_level)
}
