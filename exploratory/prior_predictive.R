# boilerplate
library(tidyverse)
library(phyloseq)
library(stray)
library(ROL)
library(driver)
library(matrixsampling)

# periodic kernel from the Kernel Cookbook
PER <- function(X, sigma = 1, rho = 1, period = 24, jitter = 0) {
  dist <- as.matrix(dist(t(X)))
  G <- sigma^2 * exp(-2*(sin(pi*dist/period)^2)/(rho^2)) + jitter*diag(ncol(dist))
  return(G)
}

build_design_matrix <- function(metadata, metadata.diet = NULL) {
  baseline_date <- metadata$collection_date[1]
  X <- sapply(metadata$collection_date, function(x) round(difftime(x, baseline_date, units="days"))) + 1
  # render X a row matrix
  dim(X) <- c(1, length(X))
  
  if(!is.null(metadata.diet)) {
    # just a few diet PCs for now
    X <- rbind(X, t(metadata.diet[,c("diet_PC1","diet_PC2","diet_PC3"),drop=F]))
    X <- rbind(X, t(metadata.diet[,c("rain_monthly","tempmax_monthly")]))
  }
  
  # do samples line up?
  if(sum((metadata$sample_id == metadata.diet$sample_id) == FALSE) > 0) {
    stop("Metadata and metadata-diet samples don't line up!\n")
  }
  
  # standardize diet and climate predictors and impute NAs
  if(nrow(X) > 1) {
    for(i in 2:nrow(X)) {
      X[i,] <- scale(X[i,])
      if(length(which(is.na(X[i,]))) > 0) {
        x <- 1:length(X[i,])
        y <- X[i,]
        X[i,] <- approx(x = x, y = y, xout = x, method = "linear")$y
      }
    }
  }
  return(X)
}

# pulled from ROL code
formalize_parameters <- function(data) {
  metadata <- sample_data(data)
  sample_status <- metadata$sample_status # replicate labels
  
  # pick an intermediate abundance ALR reference taxon
  counts <- otu_table(data)@.Data
  alr_ref <- which(order(apply(counts, 2, mean)) == round(ncol(counts)/2))
  
  # apply the ALR we'll model in
  log_ratios <- alr(counts + 0.5, d=alr_ref) # samples x taxa
  
  # get mean within-host ALR total variance
  host_var <- c()
  for(host in unique(metadata$sname)) {
    retain_sids <- metadata[metadata$sname == host,]$sample_id
    if(length(retain_sids) > 3) {
      lr <- log_ratios[(rownames(log_ratios) %in% retain_sids),]
      lr <- scale(lr, scale=FALSE)
      sample_var <- cov(t(lr))
      host_var <- c(host_var, diag(sample_var))
    }
  }
  
  mean_total_var <- mean(host_var)
  
  # return squared exponential kernel variance as 45% (90%/2) the total empirical ALR variance and the periodic
  # kernel variance as 5% (10%/2) the total empirical ALR variance
  return(list(total_variance=mean_total_var, alr_ref=alr_ref))
}

# visualize posterior intervals over Lambda, Eta
# X_predict and X_train are expected to be numeric atomic vectors (containing time points)
posterior_plot <- function(X_predict, Y_predict, taxon_idx, X_train = NULL, Y_train = NULL) {
  sample_set <- Y_predict[taxon_idx,,]
  sample_df <- gather_array(sample_set, "logratio_value", "observation", "sample_number")
  sample_df$observation <- sample_df %>%
    pull(observation) %>%
    plyr::mapvalues(., 1:nrow(sample_set), X_predict)
  quantiles <- as.data.frame(sample_df %>%
    group_by(observation) %>%
    summarise(p2.5 = quantile(logratio_value, prob=0.025),
              p25 = quantile(logratio_value, prob=0.25),
              mean = mean(logratio_value),
              p75 = quantile(logratio_value, prob=0.75),
              p97.5 = quantile(logratio_value, prob=0.975)) %>%
    ungroup())
  
  if(!is.null(X_train) & !is.null(Y_train)) {
    lr_tidy <- data.frame(observation = c(X_train), logratio_value = Y_train[taxon_idx,])
  }
  
  p <- ggplot(quantiles, aes(x=observation, y=mean)) +
    geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
    geom_ribbon(aes(ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
    geom_line(color="blue")
  if(!is.null(X_train) & !is.null(Y_train)) {
    p <- p + geom_point(data=lr_tidy, aes(x=observation, y=logratio_value), alpha=0.5)
  }
  p <- p + theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle=45)) +
    ylab("LR coord")
  return(p)
}

host <- "ACA"
cat("Fitting host:",host,"\n")

# load data
data <- load_data(tax_level="ASV")
params <- formalize_parameters(data)
data <- subset_samples(data, sname == host)
metadata <- sample_data(data)

# read diet covariate data
data.diet <- readRDS("input/ps_w_covs.RDS")
data.name_mapping <- read.csv("input/host_subject_id_to_sname_key.csv")
data.name_mapping <- unique(data.name_mapping[,c("sname","host_subject_id2")])
host.num <- as.character(data.name_mapping[data.name_mapping$sname == host,]$host_subject_id2)
data.diet <- subset_samples(data.diet, host == host.num)
metadata.diet <- sample_data(data.diet)

# pull out the count table (taxa x samples)
Y <- t(otu_table(data)@.Data)

# pull dimensions (again)
D <- nrow(Y)
N <- ncol(Y)

# ALR prior for Sigma (bacterial covariance)
upsilon <- D - 1 + 10 # specify low certainty/concentration
GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference
Xi <- GG%*%(diag(D))%*%t(GG) # take diag as covariance over log abundances
Xi <- Xi*(upsilon-D-1)

# define the prior over baselines as the empirical mean alr(Y)
alr_ys <- driver::alr((t(Y) + 0.5))
alr_means <- colMeans(alr_ys)
Theta <- function(X) matrix(alr_means, D-1, ncol(X))

X.basic <- build_design_matrix(metadata)
X.extra <- build_design_matrix(metadata, metadata.diet)

# strip off sequence variant labels
colnames(Y) <- NULL
rownames(Y) <- NULL

Y <- Y[c(setdiff(1:D,params$alr_ref),params$alr_ref),]

# define the composite kernel over samples
Gamma.basic <- function(X) {
  dc <- 0.1 # desired minimum correlation
  rho <- sqrt(-90^2/(2*log(dc))) # back calculate the decay (90 days to a drop-off of)
  Gamma_scale <- params$total_variance
  Gamma_scale <- 1
  SE(X, sigma = sqrt(Gamma_scale), rho = rho, jitter = 1e-08)
}

# define the composite kernel over samples
# FYI (myself): the optimization is pretty brittle here to the scale of Gamma giving failed (?) convergence
#               (Cholesky of Hessian fails)
Gamma.extra <- function(X) {
  jitter <- 1e-08
  # back calculate the decay in correlation to approx. 0.1 at 90 days
  dc <- 0.1 # desired minimum correlation
  rho <- sqrt(-90^2/(2*log(dc))) # back calculate the decay (90 days to a drop-off of)
  Gamma_scale <- params$total_variance
  Gamma_scale <- 1
  base_sigma_sq <- Gamma_scale * 0.5
  PER_sigma_sq <- Gamma_scale * 0.5
  SE(X[1,,drop=F], sigma = sqrt(base_sigma_sq), rho = rho, jitter = jitter) +
    PER(X[1,,drop=F], sigma = sqrt(PER_sigma_sq/2), rho = 1, period = 365, jitter = jitter) +
    SE(X[2:6,,drop=F], sigma = sqrt(PER_sigma_sq/2), rho = 1, jitter = jitter) # a few diet PCs
}

# ------------------------------------------------------------------------------------------------------
#   SIMULATE ETA FROM SE-ONLY KERNEL  -- DO THESE DATA HAVE THE APPROPRIATE VARIANCE CHARACTERISTICS?
# ------------------------------------------------------------------------------------------------------

N <- 365*2
X.test <- t(1:N)
posterior_samples <- 25
Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
Lambda <- rmatrixnormal(posterior_samples, Theta(X.test), Sigma, Gamma.basic(X.test))
Eta <- array(NA, dim = dim(Lambda))
for(i in 1:dim(Lambda)[3]) {
  Eta[,,i] <- rmatrixnormal(1, Lambda[,,i], Sigma, diag(N))[,,1]
}
# looks OK I think; the model hasn't "seen" the data at this point, so the trajectory should be
# (1) basically flat and (2) basically include the data points
taxon_idx <- 1
posterior_plot(X.test[1,], Eta, taxon_idx = taxon_idx)
# plot(X.basic[1,], t(alr_ys)[1,])

# does this have the appropriate autocorrelation structure?
posterior_sample <- sample(1:posterior_samples)[1]
centered_Eta <- t(apply(Eta[,,posterior_sample], 1, function(x) scale(x, center = TRUE, scale = FALSE)))
lags <- list()
max_lag <- 180
for(i in 1:(max_lag-1)) {
  for(j in (i+1):max_lag) {
    lag <- as.character(abs(X.test[1,i] - X.test[1,j]))
    ac <- cor(centered_Eta[,i], centered_Eta[,j])
    if(exists(lag, where = lags)) {
      lags[[lag]] <- c(lags[[lag]], ac)
    } else {
      lags[[lag]] <- c(ac)
    }
  }
}
mean_ac <- sapply(lags, mean)
plot(as.numeric(names(lags))[1:max_lag], mean_ac[1:max_lag], ylim=c(0, 1))

# ------------------------------------------------------------------------------------------------------
#   SIMULATE ETA FROM COMPOSITE KERNEL -- DO THESE DATA HAVE THE APPROPRIATE VARIANCE CHARACTERISTICS?
# ------------------------------------------------------------------------------------------------------

N <- 365*2
X.test <- t(1:N)
# simulate diet and climate components by taking real components in this time range and interpolating
selected_idx <- which(X.extra[1,] <= N)
n_features <- nrow(X.extra)-1
for(f in 1:n_features) {
  interpolated <- approx(x = X.extra[1,selected_idx], y = X.extra[f+1,selected_idx], xout = X.test[1,], method = "linear")$y
  # catch NA's at end
  last_val <- interpolated[max(which(!is.na(interpolated)))]
  interpolated[is.na(interpolated)] <- last_val
  X.test <- rbind(X.test, t(interpolated))
}

posterior_samples <- 25
Sigma <- rinvwishart(1, upsilon, Xi)[,,1]
Lambda <- rmatrixnormal(posterior_samples, Theta(X.test), Sigma, Gamma.extra(X.test))
Eta <- array(NA, dim = dim(Lambda))
for(i in 1:dim(Lambda)[3]) {
  Eta[,,i] <- rmatrixnormal(1, Lambda[,,i], Sigma, diag(N))[,,1]
}
# looks OK I think; the model hasn't "seen" the data at this point, so the trajectory should be
# (1) basically flat and (2) basically include the data points
taxon_idx <- 1
posterior_plot(X.test[1,], Eta, taxon_idx = taxon_idx)
# plot(X.basic[1,], t(alr_ys)[1,])

# does this have the appropriate autocorrelation structure?
posterior_sample <- sample(1:posterior_samples)[1]
centered_Eta <- t(apply(Eta[,,posterior_sample], 1, function(x) scale(x, center = TRUE, scale = FALSE)))
lags <- list()
max_lag <- 180
for(i in 1:(max_lag-1)) {
  for(j in (i+1):max_lag) {
    lag <- as.character(abs(X.test[1,i] - X.test[1,j]))
    ac <- cor(centered_Eta[,i], centered_Eta[,j])
    if(exists(lag, where = lags)) {
      lags[[lag]] <- c(lags[[lag]], ac)
    } else {
      lags[[lag]] <- c(ac)
    }
  }
}
mean_ac <- sapply(lags, mean)
plot(as.numeric(names(lags))[1:max_lag], mean_ac[1:max_lag], ylim=c(-0.5, 1))












