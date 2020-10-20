# simulate from GP with covariates and re-fit this with basset to judge performance

rm(list = ls())

setwd("C:/Users/kim/Documents/ROL/exploratory")

library(fido)
library(dplyr)
library(matrixsampling)
library(driver)
library(ggplot2)
library(RColorBrewer)
library(matrixsampling)

# =================================================================================================
#   VISUALIZATIONS, ETC.
# =================================================================================================

# `data` is a matrix; can reasonably be counts, log counts, Eta, or Lambda
plot_trajectory <- function(data) {
  df <- data.frame(x = c(), y = c(), taxon = c())
  for(i in 1:4) {
    df <- rbind(df, data.frame(x = 1:ncol(data), y = data[i,], taxon = i))
  }
  df$taxon <- as.factor(df$taxon)
  p <- ggplot(df) +
    geom_line(aes(x = x, y = y)) +
    facet_wrap(vars(taxon), nrow = 2, ncol = 2)
  show(p)
  return(p)
}

plot_bars <- function(counts) {
  props <- apply(counts, 2, function(x) { x/sum(x) } )
  df <- driver::gather_array(props, "proportion", "taxon", "sample")
  categories <- unique(df$taxon)
  coul = brewer.pal(4, "Spectral")
  coul = colorRampPalette(coul)(length(unique(df$taxon)))
  df$taxon <- as.factor(df$taxon)
  p <- ggplot(df, aes(x = sample, y = proportion, fill = taxon)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = coul)
  show(p)
  return(p)
}

plot_ac <- function(counts, lag_upper_limit = NULL) {
  alr.ys <- alr(t(counts) + 0.5) # samples x logratios
  alr.ys <- scale(alr.ys, center = TRUE, scale = FALSE)
  alr.ys <- t(alr.ys) # logratios x samples
  N <- ncol(alr.ys)
  lagged_error <- data.frame(lag = c(), cor = c())
  if(is.null(lag_upper_limit)) {
    lag_upper_limit <- N
  }
  for(i in 1:(lag_upper_limit-1)) {
    for(j in 2:lag_upper_limit) {
      lagged_error <- rbind(lagged_error, data.frame(lag = abs(i - j), correlation = cor(alr.ys[,i], alr.ys[,j])))
    }
  }
  plot_df <- as.data.frame(lagged_error %>%
                             group_by(lag) %>%
                             summarize(mean = mean(correlation)))
  p <- ggplot(plot_df) +
    geom_line(aes(x = lag, y = mean))
  show(p)
  return(p)
}

# periodic kernel from the Kernel Cookbook
PER <- function(X, sigma = 1, rho = 1, period = 24, jitter = 0) {
  dist <- as.matrix(dist(t(X)))
  G <- sigma^2 * exp(-2*(sin(pi*dist/period)^2)/(rho^2)) + jitter*diag(ncol(dist))
  return(G)
}

Riemann_dist <- function(A, B) {
  cholA <- t(chol(A))
  cholAInv <- solve(cholA)
  X <- cholAInv%*%B%*%(t(cholAInv))
  return(sqrt(sum(abs(log(eigen(X)$values)))))
}

Euclidean_dist <- function(A, B) {
  n <- nrow(A)*ncol(A)
  return(sqrt(sum((c(A) - c(B))^2))/n)
}

# =================================================================================================
#   SIMULATION 1: PERIODIC/OSCILLATORY TREND
# =================================================================================================

simulate.1 <- function(Sigma_scale = 1, Gamma_scale = 1, Eta_noise_scale = 1,
                       alpha = 0.5, period = 7, D = 10, cycles = 2) {
  N <- period*cycles
  X <- rbind(1:N)
  upsilon <- (D-1) + 10
  Xi <- (matrix(.4, D-1, D-1) + diag(D-1))*(upsilon - (D-1) - 1)*Sigma_scale
  Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
  baselines.true <- matrix(rnorm(D-1, 0, 1), D-1, 1)
  
  Theta.true <- function(X) {
    # simulate a treatment effect
    matrix(baselines.true, length(baselines.true), ncol(X))
  }
  
  Gamma.true <- function(X) {
    jitter <- 1e-08
    dc <- 0.1 # desired minimum correlation
    SE_days_to_baseline <- 7
    SE_rho <- sqrt(-SE_days_to_baseline^2/(2*log(dc))) # back calculate the decay
    # otherwise, simulate mixed SE + PER
    var.SE <- Gamma_scale*(1-alpha)
    var.PER <- Gamma_scale*alpha
    SE(X, sigma = sqrt(var.SE), rho = SE_rho, jitter = jitter) + PER(X, sigma = sqrt(var.PER), rho = 1, period = period, jitter = jitter)
  }
  
  Lambda <- rmatrixnormal(1, Theta.true(X), Sigma, Gamma.true(X))[,,1]
  Eta <- rmatrixnormal(1, Lambda, Sigma, diag(N)*Eta_noise_scale)[,,1]
  
  # draw sampled counts
  Y <- matrix(NA, D, N)
  for(nn in 1:N) {
    pi.nn <- alrInv(Eta[,nn])
    Y[,nn] <- rmultinom(1, rpois(1, 5000), pi.nn)
  }
  
  return(list(X = X, Y = Y, Eta = Eta, Lambda = Lambda, Sigma = Sigma,
              Theta = Theta.true, Gamma = Gamma.true, upsilon = upsilon, Xi = Xi))
}

# Y is a count matrix of taxa x samples
fit_mixed_kernels.1 <- function(data, alpha, period, Gamma_scale = 1) {
  Y <- data$Y
  N <- ncol(Y)
  X <- data$X
  upsilon <- data$upsilon
  Xi <- data$Xi
  Sigma <- data$Sigma
  # empirical estimate of baselines for each taxon
  baselines.infer <- alr(t(Y) + 0.5)
  baselines.infer <- colMeans(baselines.infer)
  Theta.infer <- function(X) {
    # no baseline effect of time or treatment
    matrix(baselines.infer, length(baselines.infer), ncol(X))
  }
  
  Gamma.infer <- function(X) {
    jitter <- 1e-08
    dc <- 0.1 # desired minimum correlation
    SE_days_to_baseline <- 7
    SE_rho <- sqrt(-SE_days_to_baseline^2/(2*log(dc))) # back calculate the decay
    # otherwise, simulate mixed SE + PER
    var.SE <- Gamma_scale*(1-alpha)
    var.PER <- Gamma_scale*alpha
    SE(X, sigma = sqrt(var.SE), rho = SE_rho, jitter = jitter) + PER(X, sigma = sqrt(var.PER), rho = 1, period = period, jitter = jitter)
  }
  fit <- fido::basset(Y, X, upsilon, Theta.infer, Gamma.infer, Xi)
  return(list(Eta = apply(fit$Eta, c(1,2), mean), Lambda = apply(fit$Lambda, c(1,2), mean),
              Sigma = apply(fit$Sigma, c(1,2), mean), Theta = Theta.infer(X), Gamma = Gamma.infer(X)))
}

# these parameter choices give a something that looks about like the empirical data
set.seed(2)
period <- 30
cycles <- 4
N <- period*cycles
data <- simulate.1(Sigma_scale = 1, Gamma_scale = 1, Eta_noise_scale = 1, alpha = 0.5, period = period, cycles = cycles)
p <- plot_ac(data$Y, lag_upper_limit = round(N*0.75))
ggsave("sim_1_autocorrelation.png", p, units = "in", dpi = 100, width = 6, height = 4)
p <- plot_bars(data$Y)
ggsave("sim_1_stackedbars.png", p, units = "in", dpi = 100, width = 9, height = 4)

# compare fits with and without simulated extra component
# fit with NO periodic component
fit.a <- fit_mixed_kernels.1(data, alpha = 0, period = period)
# fit with 50% periodic component (as simulated); tends to be a little better but it's subtle
fit.b <- fit_mixed_kernels.1(data, alpha = 0.5, period = period)

png("sim_1_Sigmas.png", height = 300, width = 850)
par(mfrow=c(1,3))
image(fit.a$Sigma)
image(data$Sigma)
image(fit.b$Sigma)
dev.off()

# model A ("bad")
Riemann_dist(data$Sigma, fit.a$Sigma)
Euclidean_dist(data$Sigma, fit.a$Sigma)

# model B ("good")
Riemann_dist(data$Sigma, fit.b$Sigma)
Euclidean_dist(data$Sigma, fit.b$Sigma)

# =================================================================================================
#   SIMULATION 2: BEFORE/AFTER "TREATMENT" EFFECT
# =================================================================================================

simulate.2 <- function(Sigma_scale = 1, Gamma_scale = 1, Eta_noise_scale = 1, binary_effect = 0, alpha = 0.5, D = 10) {
  N <- 100
  X <- rbind(1:N)
  upsilon <- (D-1) + 10
  Xi <- (matrix(.4, D-1, D-1) + diag(D-1))*(upsilon - (D-1) - 1)*Sigma_scale
  Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
  baselines.true <- matrix(rnorm(D-1, 0, 1), D-1, 1)
  
  Theta.true <- function(X) {
    # simulate a treatment effect
    matrix(baselines.true, length(baselines.true), ncol(X))
  }
  
  Gamma.true <- function(X) {
    fixed_kernel <- matrix(-binary_effect, N, N)
    fixed_kernel[1:(N/2),1:(N/2)] <- binary_effect
    fixed_kernel[(N/2+1):N,(N/2+1):N] <- binary_effect
    diag(fixed_kernel) <- 1
    jitter <- 1e-08
    dc <- 0.1 # desired minimum correlation
    SE_days_to_baseline <- 7
    SE_rho <- sqrt(-SE_days_to_baseline^2/(2*log(dc))) # back calculate the decay
    var.SE <- Gamma_scale*(1-alpha)
    var.fixed <- Gamma_scale*alpha
    SE(X[1,,drop=F], sigma = sqrt(var.SE), rho = SE_rho, jitter = jitter) + var.fixed*fixed_kernel
  }
  
  Lambda <- rmatrixnormal(1, Theta.true(X), Sigma, Gamma.true(X))[,,1]
  Eta <- rmatrixnormal(1, Lambda, Sigma, diag(N)*Eta_noise_scale)[,,1]
  
  # draw sampled counts
  Y <- matrix(NA, D, N)
  for(nn in 1:N) {
    pi.nn <- alrInv(Eta[,nn])
    Y[,nn] <- rmultinom(1, rpois(1, 5000), pi.nn)
  }
  
  return(list(X = X, Y = Y, Eta = Eta, Lambda = Lambda, Sigma = Sigma,
              Theta = Theta.true, Gamma = Gamma.true, upsilon = upsilon, Xi = Xi))
}

# Y is a count matrix of taxa x samples
fit_mixed_kernels.2 <- function(data, alpha, period, Gamma_scale = 1, binary_effect = 0) {
  Y <- data$Y
  N <- ncol(Y)
  X <- data$X
  upsilon <- data$upsilon
  Xi <- data$Xi
  Sigma <- data$Sigma
  # empirical estimate of baselines for each taxon
  baselines.infer <- alr(t(Y) + 0.5)
  baselines.infer <- colMeans(baselines.infer)
  Theta.infer <- function(X) {
    # no baseline effect of time or treatment
    matrix(baselines.infer, length(baselines.infer), ncol(X))
  }
  
  Gamma.infer <- function(X) {
    fixed_kernel <- matrix(-binary_effect, N, N)
    fixed_kernel[1:(N/2),1:(N/2)] <- binary_effect
    fixed_kernel[(N/2+1):N,(N/2+1):N] <- binary_effect
    diag(fixed_kernel) <- 1
    jitter <- 1e-08
    dc <- 0.1 # desired minimum correlation
    SE_days_to_baseline <- 7
    SE_rho <- sqrt(-SE_days_to_baseline^2/(2*log(dc))) # back calculate the decay
    var.SE <- Gamma_scale*(1-alpha)
    var.fixed <- Gamma_scale*alpha
    SE(X[1,,drop=F], sigma = sqrt(var.SE), rho = SE_rho, jitter = jitter) + var.fixed*fixed_kernel
  }
  fit <- fido::basset(Y, X, upsilon, Theta.infer, Gamma.infer, Xi)
  return(list(Eta = apply(fit$Eta, c(1,2), mean), Lambda = apply(fit$Lambda, c(1,2), mean),
              Sigma = apply(fit$Sigma, c(1,2), mean), Theta = Theta.infer(X), Gamma = Gamma.infer(X)))
}

# these parameter choices give a something that looks about like the empirical data
set.seed(1)
N <- 100 # sample size is fixed here
data <- simulate.2(Sigma_scale = 1, Gamma_scale = 1, Eta_noise_scale = 1, binary_effect = 0.8, alpha = 0.8)
p <- plot_bars(data$Y)
ggsave("sim_2_stackedbars.png", p, units = "in", dpi = 100, width = 9, height = 4)

# compare fits with and without simulated extra component
# fit with NO periodic component
fit.a <- fit_mixed_kernels.2(data, binary_effect = 0.8, alpha = 0)
# fit with 50% periodic component (as simulated); tends to be a little better but it's subtle
fit.b <- fit_mixed_kernels.2(data, binary_effect = 0.8, alpha = 0.8)

dev.off()
png("sim_2_Sigmas.png", height = 300, width = 850)
par(mfrow=c(1,3))
image(fit.a$Sigma)
image(data$Sigma)
image(fit.b$Sigma)
dev.off()

# model A ("bad")
Riemann_dist(data$Sigma, fit.a$Sigma)
Euclidean_dist(data$Sigma, fit.a$Sigma)

# model B ("good")
Riemann_dist(data$Sigma, fit.b$Sigma)
Euclidean_dist(data$Sigma, fit.b$Sigma)

# =================================================================================================
#   SIMULATION 3: SCATTERED BINARY EFFECT
# =================================================================================================

# simulate the second regime interspersed with the AC trend
simulate.3 <- function(Sigma_scale = 1, Gamma_scale = 1, Eta_noise_scale = 1, binary_effect = 0.5, alpha = 0.5, D = 10) {
  N <- 100
  X <- rbind(1:N)
  upsilon <- (D-1) + 10
  Xi <- (matrix(.4, D-1, D-1) + diag(D-1))*(upsilon - (D-1) - 1)*Sigma_scale
  Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
  baselines.true <- matrix(rnorm(D-1, 0, 1), D-1, 1)
  
  Theta.true <- function(X) {
    # simulate a treatment effect
    matrix(baselines.true, length(baselines.true), ncol(X))
  }
  
  Gamma.true <- function(X) {
    fixed_kernel <- matrix(-binary_effect, N, N)
    odds <- seq(from = 1, to = N, by = 2)
    evens <- seq(from = 2, to = N, by = 2)
    fixed_kernel[odds,odds] <- binary_effect
    fixed_kernel[evens,evens] <- binary_effect
    diag(fixed_kernel) <- 1
    jitter <- 1e-08
    dc <- 0.1 # desired minimum correlation
    SE_days_to_baseline <- 7
    SE_rho <- sqrt(-SE_days_to_baseline^2/(2*log(dc))) # back calculate the decay
    var.SE <- Gamma_scale*(1-alpha)
    var.fixed <- Gamma_scale*alpha
    SE(X[1,,drop=F], sigma = sqrt(var.SE), rho = SE_rho, jitter = jitter) + var.fixed*fixed_kernel
  }
  
  Lambda <- rmatrixnormal(1, Theta.true(X), Sigma, Gamma.true(X))[,,1]
  Eta <- rmatrixnormal(1, Lambda, Sigma, diag(N)*Eta_noise_scale)[,,1]
  
  # draw sampled counts
  Y <- matrix(NA, D, N)
  for(nn in 1:N) {
    pi.nn <- alrInv(Eta[,nn])
    Y[,nn] <- rmultinom(1, rpois(1, 5000), pi.nn)
  }
  
  return(list(X = X, Y = Y, Eta = Eta, Lambda = Lambda, Sigma = Sigma,
              Theta = Theta.true, Gamma = Gamma.true, upsilon = upsilon, Xi = Xi))
}

# Y is a count matrix of taxa x samples
fit_mixed_kernels.3 <- function(data, alpha, Gamma_scale = 1, binary_effect = 0) {
  Y <- data$Y
  N <- ncol(Y)
  X <- data$X
  upsilon <- data$upsilon
  Xi <- data$Xi
  Sigma <- data$Sigma
  # empirical estimate of baselines for each taxon
  baselines.infer <- alr(t(Y) + 0.5)
  baselines.infer <- colMeans(baselines.infer)
  Theta.infer <- function(X) {
    # no baseline effect of time or treatment
    matrix(baselines.infer, length(baselines.infer), ncol(X))
  }
  
  Gamma.infer <- function(X) {
    fixed_kernel <- matrix(-binary_effect, N, N)
    odds <- seq(from = 1, to = N, by = 2)
    evens <- seq(from = 2, to = N, by = 2)
    fixed_kernel[odds,odds] <- binary_effect
    fixed_kernel[evens,evens] <- binary_effect
    diag(fixed_kernel) <- 1
    jitter <- 1e-08
    dc <- 0.1 # desired minimum correlation
    SE_days_to_baseline <- 7
    SE_rho <- sqrt(-SE_days_to_baseline^2/(2*log(dc))) # back calculate the decay
    var.SE <- Gamma_scale*(1-alpha)
    var.fixed <- Gamma_scale*alpha
    SE(X[1,,drop=F], sigma = sqrt(var.SE), rho = SE_rho, jitter = jitter) + var.fixed*fixed_kernel
  }
  fit <- fido::basset(Y, X, upsilon, Theta.infer, Gamma.infer, Xi)
  return(list(Eta = apply(fit$Eta, c(1,2), mean), Lambda = apply(fit$Lambda, c(1,2), mean),
              Sigma = apply(fit$Sigma, c(1,2), mean), Theta = Theta.infer(X), Gamma = Gamma.infer(X)))
}

# these parameter choices give a something that looks about like the empirical data
set.seed(1)
N <- 100 # sample size is fixed here
data <- simulate.3(Sigma_scale = 1, Gamma_scale = 1, Eta_noise_scale = 1, binary_effect = 0.8, alpha = 0.8)
p <- plot_bars(data$Y)
ggsave("sim_3_stackedbars.png", p, units = "in", dpi = 100, width = 9, height = 4)

# compare fits with and without simulated extra component
# fit with NO periodic component
fit.a <- fit_mixed_kernels.3(data, binary_effect = 0.8, alpha = 0)
# fit with 50% periodic component (as simulated); tends to be a little better but it's subtle
fit.b <- fit_mixed_kernels.3(data, binary_effect = 0.8, alpha = 0.8)

dev.off()
png("sim_3_Sigmas.png", height = 300, width = 850)
par(mfrow=c(1,3))
image(fit.a$Sigma)
image(data$Sigma)
image(fit.b$Sigma)
dev.off()

# model A ("bad")
Riemann_dist(data$Sigma, fit.a$Sigma)
Euclidean_dist(data$Sigma, fit.a$Sigma)

# model B ("good")
Riemann_dist(data$Sigma, fit.b$Sigma)
Euclidean_dist(data$Sigma, fit.b$Sigma)





