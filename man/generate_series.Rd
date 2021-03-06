% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gLV_simulations.R
\name{generate_series}
\alias{generate_series}
\title{Use a generalized stochastic Lotka-Volterra model to (mechanistically) simulate count data
This code is adapted directly from Tredennick et al., 2017}
\usage{
generate_series(
  S,
  NiTh,
  AsigE,
  rm,
  al,
  K,
  sigO,
  start_time,
  end_time,
  perturbations_global,
  perturbation_sd,
  alpha
)
}
\arguments{
\item{S}{number of species to simulate}

\item{NiTh}{parameter returned by generate_innate_params (near-steady state abundances?)}

\item{AsigE}{parameter returned by generate_env_link (per-species response to environmental perturbation)}

\item{rm}{parameter returned by generate_env_link (per-species growth rates)}

\item{al}{parameter returned by generate_env_link (cooperative/competitive dynamics)}

\item{K}{per-species carrying capacities}

\item{sigO}{magnitude of observational noise}

\item{start_time}{time point at which to begin "recording" the series (previous time points will be
discarded as burn-in)}

\item{end_time}{full length of series to simulate (including burn-in)}

\item{perturbations_global}{global environmental perturbation matrix}

\item{alpha}{proportion of perturbations_global in incorporate into the local perturbations matrix}

\item{perturbations_sd}{standard deviation of zero-mean perturbations}
}
\value{
named list of simulated species series and generated parameters
}
\description{
Use a generalized stochastic Lotka-Volterra model to (mechanistically) simulate count data
This code is adapted directly from Tredennick et al., 2017
}
