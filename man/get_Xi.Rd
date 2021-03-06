% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_fitting.R
\name{get_Xi}
\alias{get_Xi}
\title{Set up a basic ALR prior}
\usage{
get_Xi(D, total_variance = 1)
}
\arguments{
\item{D}{number of features including reference (where the ALR will represent D-1)}

\item{total_variance}{scale of the log variance}
}
\value{
list containing inverse Wishart parameters degrees of freedom and scale matrix
}
\description{
Set up a basic ALR prior
}
\examples{
params <- get_Xi(D=100, log_var_scale=1)
Sigma <- matrixsampling::rinvwishart(1, params$upsilon, params$Xi)[,,1]
}
