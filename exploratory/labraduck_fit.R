library(fido)
library(phyloseq)
library(driver)
library(ggplot2)

## Here we're fitting the labraduck model to a few of the best sampled hosts in the ABRP data set.
## Usage: labraduck_fit.R host_shortname covariate_flag (opt: start time idx) (opt: end time idx)

options <- commandArgs(trailingOnly = TRUE)
if(length(options) < 1) {
  stop("Missing host argument!")
}
host <- options[1]
use_covariates <- as.logical(options[2])
begin_idx <- NULL
end_idx <- NULL
if(length(options) >= 4) {
  begin_idx <- as.numeric(options[3])
  end_idx <- as.numeric(options[4])
}

## Note: for testing, ACA with observation indices 10:30 is good...
# host <- "ACA"
# begin_idx <- 10
# end_idx <- 30
# use_covariates <- TRUE

## Local
setwd("C:/Users/kim/Documents/ROL")
#setwd("/data/mukherjeelab/roche/ROL")

## -------------------------------------------------------------------------------------------------
##  Functions
## -------------------------------------------------------------------------------------------------

get_quantiles <- function(samples_obj, subject_dates, point_estimate) {
  quantiles <- NULL
  T <- nrow(samples_obj)
  for(tt in 1:T) {
    t_quantiles <- unname(quantile(samples_obj[tt,], probs = c(0.025, 0.25, 0.75, 0.975)))
    addend <- data.frame(observation = tt,
                         q2.5 = t_quantiles[1],
                         q25 = t_quantiles[2],
                         q50 = t_quantiles[3],
                         q97.5 = t_quantiles[4]
    )
    if(is.null(quantiles)) {
      quantiles <- addend
    } else {
      quantiles <- rbind(quantiles, addend)
    }
  }
  points0 <- rowMeans(point_estimate)
  points <- data.frame(observation = subject_dates, y = points0)
  return(list(quantiles = quantiles, points = points))
}

plot_Eta <- function(fit, lr_idx) {
  subject_dates <- fit$observations
  T <- max(subject_dates)
  n_samples <- dim(fit$Eta)[3]
  D <- fit$D
  Eta_samples <- matrix(NA, T, n_samples)
  for(k in 1:n_samples) {
    EtasS_1T <- fit$Eta_DLM[,k]
    dim(EtasS_1T) <- c(D-1, T)
    Eta_samples[,k] <- EtasS_1T[lr_idx,]
  }
  plot_data <- get_quantiles(Eta_samples, subject_dates, fit$Eta[lr_idx,,])
  quantiles <- plot_data$quantiles
  points <- plot_data$points
  p <- ggplot() +
    geom_ribbon(data = quantiles, aes(x = observation, ymin = q2.5, ymax = q97.5), fill = "#dddddd") +
    geom_ribbon(data = quantiles, aes(x = observation, ymin = q25, ymax = q50), fill = "#bbbbbb") +
    geom_point(data = points, aes(x = observation, y = y)) +
    xlab("") +
    ylab("")
  p
}

plot_Theta <- function(fit, lr_idx, cov_idx = 1) {
  F <- fit$F
  Q <- nrow(F)
  F <- F[cov_idx,]
  subject_dates <- fit$observations
  T <- max(subject_dates)
  n_samples <- dim(fit$Eta)[3]
  D <- fit$D
  Theta_baseline <- numeric(n_samples)
  Theta_samples <- matrix(NA, T, n_samples)
  for(k in 1:n_samples) {
    if(cov_idx > 1) {
      ThetasS_1T <- fit$Thetas_smoothed[,k]
      dim(ThetasS_1T) <- c(Q, D-1, T)
      Theta_baseline[k] <- mean(ThetasS_1T[1,lr_idx,])
    }
    ThetasS_1T <- fit$Thetas_smoothed[,k]
    dim(ThetasS_1T) <- c(Q, D-1, T)
    ThetasS_1T <- ThetasS_1T[cov_idx,lr_idx,]
    Theta_samples[,k] <- F*ThetasS_1T
  }
  # For visualization of trajectories associated with individual covariates (not including)
  # the average captured by the intercept is added and a little bit of noise (just to ensure
  # we have an interval for geom_ribbon() to work with. This is just meant to make for quick
  # visualization of the relative effects of the individual covariates.
  Theta_samples <- Theta_samples + mean(Theta_baseline) + rnorm((D-1)*T, 0, 0.1)
  plot_data <- get_quantiles(Theta_samples, subject_dates, fit$Eta[lr_idx,,])
  quantiles <- plot_data$quantiles
  points <- plot_data$points
  p <- ggplot() +
    geom_ribbon(data = quantiles, aes(x = observation, ymin = q2.5, ymax = q97.5), fill = "#dddddd") +
    geom_ribbon(data = quantiles, aes(x = observation, ymin = q25, ymax = q50), fill = "#bbbbbb") +
    geom_point(data = points, aes(x = observation, y = y)) +
    xlab("") +
    ylab("")
  p
}

pull_data <- function(host) {
  host <<- host # weird hack for phyloseq subset_samples()
  data <- readRDS("input/filtered_family_5_20.rds")
  data <- subset_samples(data, sname == host)
  md <- sample_data(data)
  
  # pull diet and climate data
  data.diet <- readRDS("input/ps_w_covs.RDS")
  data.name_mapping <- read.csv("input/host_subject_id_to_sname_key.csv")
  data.name_mapping <- unique(data.name_mapping[,c("sname","host_subject_id2")])
  host.num <<- as.character(data.name_mapping[data.name_mapping$sname == host,]$host_subject_id2)
  data.diet <- subset_samples(data.diet, host == host.num)
  metadata.diet <- sample_data(data.diet)
  
  days <- md$collection_date
  day0 <- min(days)
  days_baseline <- round(unname(sapply(days, function(x) difftime(x, day0, units = "days")))) + 1
  season <- md$season
  season[season == "Dry"] <- 0
  season[season == "Wet"] <- 1
  season <- as.numeric(season)
  # order is verifiably the same here
  diet_PC1 <- as.vector(scale(metadata.diet$PC1))
  rain_monthly <- as.vector(scale(metadata.diet$rain_monthly))
  counts <- t(as.matrix(otu_table(data)@.Data))
  rownames(counts) <- NULL
  colnames(counts) <- NULL
  return(list(counts = counts, days = days_baseline,
              season = season, diet = diet_PC1, rain = rain_monthly))
}

slice_dataset <- function(data, begin_idx, end_idx) {
  Y <- data$counts
  D <- nrow(Y)
  N <- ncol(Y)
  if(is.null(begin_idx) | is.null(end_idx)) {
    begin_idx <- 1
    end_idx <- N
  }
  # subset all
  data$counts <- Y[,begin_idx:end_idx]
  data$season <- data$season[begin_idx:end_idx]
  data$diet <- data$diet[begin_idx:end_idx]
  data$rain <- data$rain[begin_idx:end_idx]
  data$days <- data$days[begin_idx:end_idx]
  # shift days to start at 1
  data$days <- data$days - min(data$days) + 1
  return(data)
}

fit_model <- function(data, use_covariates) {
  T <- max(data$days)
  # Build the pseudo-covariate matrix
  if(use_covariates) {
    F <- matrix(0, 3, T)
    F[1,] <- 1
    for(i in 1:length(data$days)) {
      # F[2,data$days[i]] <- data$season[i] # season actually worsens the fit here
      F[2,data$days[i]] <- data$diet[i]
      F[3,data$days[i]] <- data$rain[i]
    }
  } else {
    F <- matrix(1, 1, T)
  }
  Q <- nrow(F)
  D <- nrow(data$counts)
  W <- diag(Q)
  G <- diag(Q)
  
  Y <- data$counts

  upsilon <- D-1+10
  GG <- cbind(diag(D-1), -1)
  Xi <- GG%*%(diag(D)*1)%*%t(GG)
  Xi <- Xi*(upsilon-D-1)
  
  C0 <- W
  alr_ys <- driver::alr((t(Y) + 0.5))
  alr_means <- colMeans(alr_ys)
  M0 <- matrix(0, Q, D-1)
  M0[1,] <- alr_means
  
  start <- Sys.time()
  fit <- labraduck(Y = Y, upsilon = upsilon, Xi = Xi, F = F, G = G, W = W, M0 = M0, C0 = C0,
                   observations = data$days, gamma_scale = 1, W_scale = 1, apply_smoother=TRUE)
  end <- Sys.time()
  return(list(fit = fit, runtime = end-start))
}

## -------------------------------------------------------------------------------------------------
##  Main
## -------------------------------------------------------------------------------------------------

data <- pull_data(host)
data <- slice_dataset(data, begin_idx, end_idx)

tag <- "cov"
if(!use_covariates) {
  tag <- paste0("no",tag)
}
base_fn <- paste0("fit_",tag,"_",host)

fit <- fit_model(data, use_covariates = use_covariates)
cat("Fit time:",fit$runtime,"sec\n")
fit <- fit$fit
saveRDS(fit, paste0(base_fn,".rds"))

p <- plot_Eta(fit, 1)
ggsave(paste0(base_fn,"_Theta.png"), p, dpi = 100, units = "in", height = 4, width = 12)

if(use_covariates) {
  for(covariate in 1:3) {
    p <- plot_Theta(fit, 1, covariate)
    ggsave(paste0(base_fn,"_Eta",covariate,".png"), p, dpi = 100, units = "in", height = 4, width = 10)
  }
}

