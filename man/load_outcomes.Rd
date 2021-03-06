% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_analyses.R
\name{load_outcomes}
\alias{load_outcomes}
\title{Loads fitness annotations associated with hosts from file}
\usage{
load_outcomes(hosts)
}
\arguments{
\item{hosts}{a vector of hosts by sname}
}
\value{
data.frame
}
\description{
Loads fitness annotations associated with hosts from file
}
\details{
Output is a data.frame of fitness annotations sorted by host short name
}
\examples{
outcomes <- load_outcomes(data)
}
