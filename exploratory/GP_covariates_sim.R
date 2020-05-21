# simulate from GP with covariates and re-fit this with basset to judge performance

rm(list = ls())

library(stray)
library(matrixsampling)
library(driver)
library(ggplot2)

# periodic kernel from the Kernel Cookbook
PER <- function(X, sigma = 1, rho = 1, period = 24, jitter = 0) {
  dist <- as.matrix(dist(t(X)))
  G <- sigma^2 * exp(-2*(sin(pi*dist/period)^2)/(rho^2)) + jitter*diag(ncol(dist))
  return(G)
}

# =================================================================================================
#   SIMULATE FROM MIXED PERIODIC GP
# =================================================================================================

simulate <- function(alpha = 0.5, period = 24, D = 20, N = 24*2) {
  X <- rbind(1:N)
  upsilon <- D + 10
  Xi <- (matrix(.4, D-1, D-1) + diag(D-1))*(upsilon - D - 1)
  Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
  baselines.true <- matrix(rnorm(D-1, 0, 2), D-1, 1)
  
  Theta.true <- function(X) {
    # simulate a treatment effect
    matrix(baselines.true, length(baselines.true), ncol(X))
  }
  
  Gamma.true <- function(X) {
    jitter <- 1e-08
    dc <- 0.1 # desired minimum correlation
    SE_days_to_baseline <- 7
    SE_rho <- sqrt(-SE_days_to_baseline^2/(2*log(dc))) # back calculate the decay
    total_variance <- 1
    alpha <- alpha
    var.SE <- total_variance*(1-alpha)
    var.PER <- total_variance*alpha
    SE(X, sigma = sqrt(var.SE), rho = SE_rho, jitter = jitter) + PER(X, sigma = sqrt(var.PER), rho = 1, period = 24, jitter = jitter)
  }
  
  Lambda <- rmatrixnormal(1, Theta.true(X), Sigma, Gamma.true(X))[,,1]
  Eta <- rmatrixnormal(1, Lambda, Sigma, diag(N))[,,1]
  
  # draw sampled counts
  Y <- matrix(NA, D, N)
  for(nn in 1:N) {
    pi.nn <- alrInv(Eta[,nn])
    Y[,nn] <- rmultinom(1, rpois(1, 5000), pi.nn)
  }
  
  return(list(X = X, Y = Y, Sigma = Sigma, upsilon = upsilon, Xi = Xi))
}

# =================================================================================================
#   PARAM SWEEP WITH KERNELS
# =================================================================================================

# Y is a count matrix of taxa x samples
fit_mixed_kernels <- function(data, alpha, return_model = FALSE) {
  Y <- data$Y
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
    total_variance <- 1
    var.SE <- total_variance*(1-alpha)
    var.PER <- total_variance*alpha
    SE(X, sigma = sqrt(var.SE), rho = SE_rho, jitter = jitter) + PER(X, sigma = sqrt(var.PER), rho = 1, period = 24, jitter = jitter)
  }
  fit <- stray::basset(Y, X, upsilon, Theta.infer, Gamma.infer, Xi)
  ret_obj <- list(error = sqrt(sum((c(Sigma) - c(apply(fit$Sigma, c(1,2), mean)))^2)))
  if(return_model) {
    ret_obj[["model"]] <- fit
  }
  return(ret_obj)
}

error_df <- data.frame(true_alpha = c(), test_alpha = c(), rmse = c())
true_alphas <- seq(from = 0, to = 1, by = 0.5)
true_alphas <- rep(true_alphas, 3)
for(ta in true_alphas) {
  cat("True alpha:",ta,"\n")
  data <- simulate(ta)
  test_alpha <- seq(from = 0, to = 1, by = 0.1)
  for(a in test_alpha) {
    error <- fit_mixed_kernels(data, a)
    error_df <- rbind(error_df, data.frame(true_alpha = ta, test_alpha = a, rmse = error))
  }
}

df <- error_df
df$true_alpha <- as.factor(df$true_alpha)
ggplot(df) +
  geom_smooth(aes(x = test_alpha, y = rmse, color = true_alpha), method = "loess") +
  geom_smooth(aes(x = test_alpha, y = rmse, color = true_alpha), method = "loess")

# not clear that there's any effect on accuracy of minimizing the difference between simulated
# and fitted kernels

# do the estimates of Sigma look good (to the eye) regardless?

# always seems worse to estimate an over-large periodic component
# but for subtle periodicity, there's not much advantage to including a periodic component
data <- simulate(0.2)
fit <- fit_mixed_kernels(data, 0.0, return_model = TRUE)
cat("RMSE:",fit$error,"\n")
par(mfrow=c(1,2))
image(data$Sigma)
image(apply(fit$model$Sigma, c(1,2), mean))

# =================================================================================================
#   HOW MUCH PERIODICITY MUST BE IN THE UNDERLYING SYSTEM FOR IT TO BE PRESENT IN AN AC PLOT?
# =================================================================================================

calc_ac <- function(Y, lagged_error = NULL) {
  alr.ys <- alr(t(Y) + 0.5) # samples x logratios
  alr.ys <- scale(alr.ys, center = TRUE, scale = FALSE)
  alr.ys <- t(alr.ys) # logratios x samples
  N <- ncol(alr.ys)
  if(is.null(lagged_error)) {
    lagged_error <- data.frame(lag = c(), cor = c())
  }
  for(i in 1:(N-1)) {
    for(j in 2:N) {
      lagged_error <- rbind(lagged_error, data.frame(lag = abs(i - j), correlation = cor(alr.ys[,i], alr.ys[,j])))
    }
  }
  return(lagged_error)
}

alpha <- 0.5
data_sets <- list()
lagged_error <- NULL
for(dd in 1:20) {
  cat("Rep:",dd,"\n")
  data <- simulate(alpha)
  lagged_error <- calc_ac(data$Y, lagged_error)
}
plot_df <- as.data.frame(lagged_error %>%
                           group_by(lag) %>%
                           summarize(mean = mean(correlation)))
ggplot(plot_df) +
  geom_point(aes(x = lag, y = mean))











