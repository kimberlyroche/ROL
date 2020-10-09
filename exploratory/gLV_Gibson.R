# This simulates absolute abundances out of a generalized Lotka-Volterra model.

library(driver)
library(ggplot2)
library(matrixsampling)
library(mvtnorm)
library(stray)
library(grid)
library(gridExtra)
library(ROL)
library(Rcpp)
library(vars)

# Yanked from universality_score_heatmap.R
calc_xy <- function(vSigmas) {
  if(is.null(dim(vSigmas))) {
    dim(vSigmas) <- c(length(vSigmas), 1)
  }
  shared_direction <- mean(apply(sign(vSigmas), 2, function(z) {
    max(table(z)/length(z))
  }))
  mean_strength <- mean(abs(apply(vSigmas, 2, mean)))
  list(x = shared_direction, y = mean_strength)
}

detrend_series <- function(X) {
  detrended_series <- X
  for(i in 1:nrow(detrended_series)) {
    detrended_series[i,] <- arima(detrended_series[i,], order = c(1, 0, 0))$residuals
  }
  return(detrended_series)
}

simulate_gLV <- function(L, T, noise_scale = 1, a1 = NULL, a2 = NULL, B = NULL, R = NULL, correlated_noise = FALSE) {
  # Initialize the count matrix and the scale (this is the starting abundance of the average feature)
  abundance_scale <- 100
  X <- matrix(0, L, T)
  X[,1] <- rnorm(L, abundance_scale, 10)

  # Growth and carrying capacity parameters
  if(is.null(a1)) {
    a1 <- runif(L, min = 0.01, max = 0.05) # growth trend
  }
  
  if(is.null(a2)) {
    a2 <- -runif(L, min = 0.01, max = 0.05)/abundance_scale # negative pressure equiv. to a carrying capacity
  }
  
  if(is.null(B)) {
    # Covariance with other taxa
    Omega <- diag(L)
    Omega <- Omega / (100 * L * abundance_scale)
    #Omega <- Omega / (100 * abundance_scale)
    nu <- 2
    B <- rinvwishart(1, L + nu, Omega*nu, checkSymmetry = FALSE)[,,1]
  }
  
  if(is.null(R)) {
    if(correlated_noise) {
      Omega <- diag(L)
      Omega <- Omega * noise_scale
      nu <- 2
      R <- rinvwishart(1, L + nu, Omega*nu, checkSymmetry = FALSE)[,,1]
    } else {
      R <- diag(L) * noise_scale
    }
  }
  
  # Simulate
  for(i in 2:T) {
    noise <- as.vector(rmvnorm(1, mean = rep(0, L), sigma = R))
    for(ll in 1:L) {
      p1 <- a1[ll]
      p2 <- a2[ll] * X[ll,i-1]
      p3 <- 0
      for(jj in setdiff(1:L, ll)) {
        p3 <- p3 + B[ll,jj] * X[jj,i-1]
      }
      p4 <- noise[ll]
      X[ll,i] <- X[ll,i-1] * (1 + p1 + p2 + p3) + p4
      if(X[ll,i] < 1) {
        # Set a floor for abundance
        X[ll,i] <- rpois(1, lambda = 1)
      }
    }
  }
  
  return(list(X = X, a1 = a1, a2 = a2, B = B, R = R))
}

visualize_covariance <- function(B) {
  df <- gather_array(cov2cor(B), "value", "taxon1", "taxon2")
  p <- ggplot(df) +
    geom_tile(aes(x = taxon1, y = taxon2, fill = value)) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)
  return(p)
}

estimated_inferred_covariance <- function(X) {
  # sample_no <- round(ncol(X)/10)
  # X.downsampled <- X[,sort(sample(1:ncol(X))[1:100])]
  # lag2_ac <- acf(X.downsampled[1,], plot = FALSE)$acf[2]
  # observed_correlation <- cov2cor(cov(t(X.downsampled)))
  observed_correlation <- cov2cor(cov(t(X)))
  # return(list(observed_correlation = observed_correlation, lag2_ac = lag2_ac))
  return(observed_correlation)
}

visualize_trajectories <- function(X) {
  df <- gather_array(X, "abundance", "taxon", "timepoint")
  df$taxon <- as.factor(df$taxon)
  p <- ggplot(df) +
    geom_path(aes(x = timepoint, y = abundance, color = taxon)) +
    geom_point(aes(x = timepoint, y = abundance, color = taxon))
  return(p)
}

visualize_proportions <- function(X) {
  p <- apply(X, 2, function(x) x/sum(x))
  df <- gather_array(p, "relative_abundance", "taxon", "timepoint")
  df$taxon <- as.factor(df$taxon)
  p <- ggplot(df, aes(fill = taxon, y = relative_abundance, x = timepoint)) +
    geom_bar(position = "stack", stat = "identity")
  return(p)
}

# Simulation 1 -- Balanced and uninteresting (little noise)
sim1 <- simulate_gLV(L = 10, T = 200, noise_scale = 1)
visualize_trajectories(sim1$X)
visualize_proportions(sim1$X)

# Simulation 2 -- Noise
sim1 <- simulate_gLV(L = 10, T = 200, noise_scale = 100)
visualize_trajectories(sim1$X)
visualize_proportions(sim1$X)

# Simulation 3 -- Correlated noise (doesn't look particularly different)
sim1 <- simulate_gLV(L = 10, T = 200, noise_scale = 100, correlated_noise = TRUE)
visualize_trajectories(sim1$X)
visualize_proportions(sim1$X)

sim1 <- simulate_gLV(L = 10, T = 200, noise_scale = 1)
sim2 <- simulate_gLV(L = 10, T = 200, noise_scale = 1, a1 = sim1$a1, a2 = sim1$a2, B = sim1$B, R = sim1$R)
grid.arrange(visualize_trajectories(sim1$X), visualize_trajectories(sim2$X), nrow = 2)

est_Sigma1 <- estimated_inferred_covariance(sim1$X[,101:200])
est_Sigma2 <- estimated_inferred_covariance(sim2$X[,101:200])

grid.arrange(visualize_covariance(sim1$B),
             visualize_covariance(sim1$R),
             visualize_covariance(est_Sigma1),
             visualize_covariance(est_Sigma2), nrow = 2, ncol = 2)

est_Sigma1 <- estimated_inferred_covariance(detrend_series(sim1$X))
est_Sigma2 <- estimated_inferred_covariance(detrend_series(sim2$X))
grid.arrange(visualize_covariance(sim1$B),
             visualize_covariance(sim1$R),
             visualize_covariance(est_Sigma1),
             visualize_covariance(est_Sigma2), nrow = 2, ncol = 2)

# Crude universality test

H <- 20
L <- 20
host_associations <- matrix(NA, H, (L^2)/2 - L/2)
a1 <- NULL
a2 <- NULL
B <- NULL
R <- NULL
for(h in 1:H) {
  cat("Simulating host",h,"\n")
  sim <- simulate_gLV(L = L, T = 200, noise_scale = 1, a1 = a1, a2 = a2, B = B, R = R)
  #sim <- simulate_gLV(L = L, T = 200, noise_scale = 100, a1 = a1, a2 = a2, B = B, R = R, correlated_noise = FALSE)
  #sim <- simulate_gLV(L = L, T = 200, noise_scale = 100, a1 = a1, a2 = a2, B = B, R = R, correlated_noise = TRUE)
  a1 <- sim$a1
  a2 <- sim$a2
  B <- sim$B
  R <- sim$R
  #temp <- estimated_inferred_covariance(sim$X[,101:200])
  temp <- estimated_inferred_covariance(sim$X[,1:100])
  host_associations[h,] <- temp[upper.tri(temp, diag = FALSE)]
}

point_data <- data.frame(x = c(), y = c())
for(i in 1:ncol(host_associations)) {
  temp <- calc_xy(host_associations[,i])
  point_data <- rbind(point_data, data.frame(x = temp$x, y = temp$y))
}

ggplot(point_data) +
  geom_point(aes(x = x, y = y)) +
  xlim(0.5, 1) +
  ylim(0, 1)

# Takeaway:
# (1) Simulated covariance in B is consistently estimable across hosts only if there's low noise and we
#     restrict ourself to the recovery-to-steady-state portion of the simulated series
# (2) Simulated covariance in R is consistently estimable across hosts once it's sufficiently large
#     relative to other effects (this represents shared reaction to perturbation, e.g. diet)




