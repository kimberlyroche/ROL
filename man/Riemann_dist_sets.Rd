% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{Riemann_dist_sets}
\alias{Riemann_dist_sets}
\title{Calculate distances between all pairs of samples from two sample sets}
\usage{
Riemann_dist_sets(A, B, n_hosts, host_samples_A, host_samples_B)
}
\arguments{
\item{A}{posterior sample matrix 1, ordered by sample number, the host concatenated across columns}

\item{B}{posterior sample matrix 2, ordered by sample number, the host concatenated across columns}

\item{n_hosts}{number of hosts in set A and B (must be equal)}

\item{host_samples_A}{number of samples per host in set A}

\item{host_samples_B}{number of samples per host in set B}
}
\value{
distance
}
\description{
Calculate distances between all pairs of samples from two sample sets
}
\details{
Dimensions of the samples matrix are D rows x DN columns where D = no. of features
and N = no. of samples. For two hosts "DUD" and "OMO" with 3 samples each, these are arranged
column-wise in the matrix as (DUD 1, OMO 1, DUD 2, OMO 2, DUD 3, OMO 3).
}
