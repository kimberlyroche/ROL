library(tidyverse)
library(phyloseq)
library(stray)
library(ROL)
#library(pracma) # for pchip()
library(driver)

args <- commandArgs(trailing = TRUE)
host <- args[1]

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
posterior_plot <- function(X_predict, Y_predict, taxon_idx, X_train = NULL, Y_train = NULL) {
  sample_set <- Y_predict[taxon_idx,,]
  sample_df <- gather_array(sample_set, "logratio_value", "observation", "sample_number")
  quantiles <- sample_df %>%
    group_by(observation) %>%
    summarise(p2.5 = quantile(logratio_value, prob=0.025),
              p25 = quantile(logratio_value, prob=0.25),
              mean = mean(logratio_value),
              p75 = quantile(logratio_value, prob=0.75),
              p97.5 = quantile(logratio_value, prob=0.975)) %>%
    ungroup()
  
  if(!is.null(X_train) & !is.null(Y_train)) {
    lr_tidy <- data.frame(observation = c(X_train[1,]), logratio_value = Y_train[taxon_idx,])
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
  Gamma_scale <- 4
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
  Gamma_scale <- 4
  base_sigma_sq <- Gamma_scale * 0.5
  PER_sigma_sq <- Gamma_scale * 0.5
  SE(X[1,,drop=F], sigma = sqrt(base_sigma_sq), rho = rho, jitter = jitter) +
    PER(X[1,,drop=F], sigma = sqrt(PER_sigma_sq/2), rho = 1, period = 365, jitter = jitter) +
    SE(X[2:6,,drop=F], sigma = sqrt(PER_sigma_sq/2), rho = 1, jitter = jitter) # a few diet PCs
}

# ALR prior for Sigma (bacterial covariance)
upsilon <- D - 1 + 10 # specify low certainty/concentration
GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference
Xi <- GG%*%(diag(D))%*%t(GG) # take diag as covariance over log abundances
Xi <- Xi*(upsilon-D-1)

# define the prior over baselines as the empirical mean alr(Y)
alr_ys <- driver::alr((t(Y) + 0.5))
alr_means <- colMeans(alr_ys)
Theta <- function(X) matrix(alr_means, D-1, ncol(X))

#cat("\tFitting basic model...\n")
# MAP fits
# fit.basic <- stray::basset(Y, X.basic, upsilon, Theta, Gamma.basic, Xi, n_samples = 0, ret_mean = TRUE)
# fit.extra <- stray::basset(Y, X.extra, upsilon, Theta, Gamma.extra, Xi, n_samples = 0, ret_mean = TRUE)
#fit.basic <- stray::basset(Y, X.basic, upsilon, Theta, Gamma.basic, Xi, n_samples = 100,
#                           b2 = 0.99, step_size = 0.004, eps_f = 1e-11, eps_g = 1e-04,
#                           max_iter = 10000L, decomp_method = "cholesky", optim_method = "lbfgs")
# saveRDS(fit.basic, file = paste0(host,"_fit_basic.rds"))

cat("\tFitting extra model...\n")
# this generally fails due to Hessian inversion problem; indicates its not finding the optimum?
# Gamma.extra(X.extra) is PSD so that's not an issue
# predictive fits look weird when it does work?
fit.extra <- stray::basset(Y, X.extra, upsilon, Theta, Gamma.extra, Xi, n_samples = 100,
                           b2 = 0.99, step_size = 0.003, eps_f = 1e-10, eps_g = 1e-04,
                           max_iter = 10000L, decomp_method = "eigen", optim_method = "lbfgs")
# saveRDS(fit.extra, file = paste0(host,"_fit_extra.rds"))

# assess comparative fit
# (1) log marginal likelihood
# write(paste0(host,"\t",fit.basic$logMarginalLikelihood,"\t",fit.extra$logMarginalLikelihood),
#             "mll_out.txt", append = TRUE)

cat("Rendering visualizations...\n")
# (2) eyeball the predictive fits for the first taxon
# basic model
X_predict.basic <- t(1:(max(X.basic)))
predicted.basic <- predict(fit.basic, X_predict.basic, response = "Eta", iter = fit.basic$iter)
plot.basic <- posterior_plot(X_predict.basic, predicted.basic, taxon_idx = 1,
              X_train = X.basic, Y_train = t(alr_ys))
# ggsave(paste0(host,"_predictive_01.png"), plot.basic, units = "in", dpi = 100, height = 4, width = 12)

# need to interpolate additional features (diet, climate) for the extra model
full_n <- max(X.extra[1,,drop=F])
partial_n <- full_n
X_predict.extra <- t(1:partial_n)
trunc_X <- X.extra[1,]
selected_idx <- trunc_X <= partial_n
trunc_X <- trunc_X[selected_idx]
n_features <- nrow(X.extra)-1
for(f in 1:n_features) {
  #X_predict.extra <- rbind(X_predict.extra, pchip(xi = trunc_X, yi = X.extra[f+1,selected_idx], x = 1:partial_n))
  X_predict.extra <- rbind(X_predict.extra, t(approx(x = trunc_X, y = X.extra[f+1,selected_idx], xout = 1:partial_n, method = "linear")$y))
}
# then predict
predicted.extra <- predict(fit.extra, X_predict.extra, response="Eta", iter=fit.extra$iter)
plot.extra <- posterior_plot(X_predict.extra, predicted.extra, taxon_idx = 1,
              X_train = X.extra[,selected_idx], Y_train = t(alr_ys)[,selected_idx])
# ggsave(paste0(host,"_predictive_02.png"), plot.extra, units = "in", dpi = 100, height = 4, width = 12)

# (3) eyeball the differences in Sigma
# png(paste0(host,"_Sigma_01.png"))
# crude_Sigma_basic <- apply(fit.basic$Sigma, c(1,2), mean)
# image(crude_Sigma_basic)
# dev.off()
# png(paste0(host,"_Sigma_02.png"))
# crude_Sigma_extra <- apply(fit.extra$Sigma, c(1,2), mean)
# image(crude_Sigma_extra)
# dev.off()

cat("Complete.\n")
