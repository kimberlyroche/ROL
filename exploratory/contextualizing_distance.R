library(ROL)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(abind)
library(matrixsampling)
library(stray)
library(driver)

## --------------------------------------------------------------------------------------------------------
##     Parse ABRP data
## --------------------------------------------------------------------------------------------------------

Sigmas <- load_MAP_estimates(tax_level = "ASV", logratio = "clr")
# convert these to correlations
Sigmas <- lapply(Sigmas, function(x) cov2cor(x)) # JOHANNES, THANK YOU
n_hosts <- length(Sigmas)
P <- dim(Sigmas[[1]])[1]
pairwise_combos <- P^2 - P^2/2 - P/2

# strip out the triangular portion of the correlation matrices, vectorize these, and put them in a matrix
# for later convenience
vectorized_Sigmas <- matrix(NA, n_hosts, pairwise_combos)
for(i in 1:n_hosts){
  vectorized_Sigmas[i,] <- Sigmas[[i]][upper.tri(Sigmas[[i]], diag = F)]
}

# get host primary social group assignments
data <- load_data(tax_level = "ASV")
labels <- list()
for(host in names(Sigmas)) {
  subdata <- subset_samples(data, sname == host)
  metadata <- sample_data(subdata)
  labels[[host]] <- names(which(table(metadata[["grp"]]) == max(table(metadata[["grp"]]))))[1]
}

## --------------------------------------------------------------------------------------------------------
##     Parse unfiltered Johnson et al. (2019) data
## --------------------------------------------------------------------------------------------------------

# read count table
counts <- read.table("input/johnson2019/taxonomy_counts_s.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# read sample ID to subject ID mapping
mapping <- read.table("input/johnson2019/SampleID_map.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

vectorized_Sigmas_johnson2019 <- NULL

for(subject in unique(mapping$UserName)) {
  cat("Evaluating subject:",subject,"\n")
  subject_samples <- mapping[mapping$UserName == subject,]$SampleID
  subject_days <- mapping[mapping$UserName == subject,]$StudyDayNo
  samples_to_pull <- colnames(counts)[which(colnames(counts) %in% subject_samples)]
  if(length(samples_to_pull) > 0) {
    subject_days <- which(subject_samples %in% samples_to_pull)
    subject_counts <- counts[,samples_to_pull] # omitting taxonomy
    dim(subject_counts)
    # later: we probably want to this about removing taxa that are very rare within any individual
    
    # fit this with stray::basset
    Y <- as.matrix(subject_counts)
    D <- nrow(Y)
    X <- matrix(subject_days, 1, length(subject_days))
    
    alr_ys <- driver::alr((t(Y) + 0.5))
    alr_means <- colMeans(alr_ys)
    Theta <- function(X) matrix(alr_means, D-1, ncol(X))
    
    taxa_covariance <- get_Xi(D, total_variance = 1)
    
    rho <- calc_se_decay(min_correlation = 0.1, days_to_baseline = 1)
    Gamma <- function(X) {
      SE(X, sigma = 1, rho = rho, jitter = 1e-08)
    }
    
    fit <- stray::basset(Y, X, taxa_covariance$upsilon, Theta, Gamma, taxa_covariance$Xi,
                         n_samples = 0, ret_mean = TRUE,
                         b2 = 0.98, step_size = 0.004, eps_f = 1e-11, eps_g = 1e-05,
                         max_iter = 10000L, optim_method = "adam")
    fit.clr <- to_clr(fit)
    fit.clr$Sigma <- fit.clr$Sigma[,,1]
    fit.clr$Sigma <- cov2cor(fit.clr$Sigma)
    if(is.null(vectorized_Sigmas_johnson2019)) {
      vectorized_Sigmas_johnson2019 <- fit.clr$Sigma[upper.tri(fit.clr$Sigma, diag = FALSE)]
    } else {
      vectorized_Sigmas_johnson2019 <- rbind(vectorized_Sigmas_johnson2019, fit.clr$Sigma[upper.tri(fit.clr$Sigma, diag = FALSE)])
    }
  }
}

# # visualize the last fitted model for plausibility
# fit <- stray::basset(Y, X, taxa_covariance$upsilon, Theta, Gamma, taxa_covariance$Xi,
#                      n_samples = 500, ret_mean = FALSE, 
#                      b2 = 0.98, step_size = 0.004, eps_f = 1e-11, eps_g = 1e-05,
#                      max_iter = 10000L, optim_method = "adam")
# Eta <- predict(fit, t(subject_days), response = "Eta", iter = fit$iter)
# Eta <- alrInv_array(Eta, fit$D, 1)
# Eta <- clr_array(Eta, 1)
# lr_ys <- clr(t(fit$Y) + 0.5)
# 
# coord <- sample(1:D)[1]
# observations <- subject_days
# lr_tidy <- gather_array(lr_ys, "logratio_value", "timepoint", "logratio_coord")
# 
# # replace timepoints with observation dates for readability
# # map <- data.frame(timepoint=1:length(observations), observation=c(observations))
# # lr_tidy <- merge(lr_tidy, map, by="timepoint")
# # lr_tidy <- lr_tidy[,!(names(lr_tidy) %in% c("timepoint"))]
# 
# no_samples <- dim(Eta)[3]
# posterior_samples <- gather_array(Eta[coord,,], "logratio_value", "observation", "sample_number")
#   
# # get quantiles
# post_quantiles <- posterior_samples %>%
#   group_by(observation) %>%
#   summarise(p2.5 = quantile(logratio_value, prob=0.025),
#             p5 = quantile(logratio_value, prob=0.05),
#             p10 = quantile(logratio_value, prob=0.1),
#             p25 = quantile(logratio_value, prob=0.25),
#             p50 = quantile(logratio_value, prob=0.5),
#             mean = mean(logratio_value),
#             p75 = quantile(logratio_value, prob=0.75),
#             p90 = quantile(logratio_value, prob=0.9),
#             p95 = quantile(logratio_value, prob=0.95),
#             p97.5 = quantile(logratio_value, prob=0.975)) %>%
#   ungroup()
#   
# p <- ggplot(post_quantiles, aes(x=observation, y=mean)) +
#   geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
#   geom_ribbon(aes(ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
#   geom_line(color="blue") +
#   geom_point(data=lr_tidy[lr_tidy$logratio_coord == coord,], aes(x=timepoint, y=logratio_value), alpha=0.5) +
#   theme_minimal() +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(angle=45)) +
#   ylab("LR coord")
# show(p)

## --------------------------------------------------------------------------------------------------------
##     Visualize ABRP social group vs. Johnson et al. (2019) human data
## --------------------------------------------------------------------------------------------------------

calc_map_xy <- function(vectorized_Sigmas) {
  within_score <- mean(apply(abs(vectorized_Sigmas), 1, mean))
  between_score <- mean(apply(combn(1:nrow(vectorized_Sigmas), m = 2), 2, function(x) {
    cor(vectorized_Sigmas[x[1],], vectorized_Sigmas[x[2],])
  }))
  return(list(x = within_score, y = between_score))
}

point <- calc_map_xy(vectorized_Sigmas_johnson2019)
plot_df <- data.frame(x = point$x, y = point$y, group = "johnson2019")
for(group in unlist(unique(labels))) {
  group_Sigmas <- vectorized_Sigmas[unname(which(labels == group)),]
  point <- calc_map_xy(group_Sigmas)
  plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = group))
}

p_combo <- ggplot(plot_df) +
  geom_point(aes(x = x, y = y, color = group), size = 3) +
  xlim(c(0,1)) +
  ylim(c(-0.1,1)) +
  xlab("avg. strength of associations (within hosts)") +
  ylab("avg. agreement between hosts")
print(p_combo)

## --------------------------------------------------------------------------------------------------------
##     Visualize on "map" as quantiles
## --------------------------------------------------------------------------------------------------------

vectorized_DUI <- vectorized_Sigmas[which(names(Sigmas) == "DUI"),]
ecdf <- sort(abs(vectorized_DUI))
# separate into quantiles
boundaries <- quantile(ecdf, probs = seq(from = 0, to = 1, by = (1/4)))

# bin a set of pairwise relationships between microbes
get_col_idx <- function(vec, lower_boundary, upper_boundary) {
  which(abs(vec) >= lower_boundary & abs(vec) < upper_boundary)
}

# col_idx are the column indices of pairs of microbes to consider here
calc_map_xy_binned <- function(vectorized_Sigmas, col_idx) {
  within_score <- mean(apply(abs(vectorized_Sigmas[,col_idx]), 1, mean))
  between_score <- mean(apply(combn(1:nrow(vectorized_Sigmas), m = 2), 2, function(x) {
    cor(vectorized_Sigmas[x[1],col_idx], vectorized_Sigmas[x[2],col_idx])
  }))
  return(list(x = within_score, y = between_score))
}

# add in baboon quantiles
plot_df <- data.frame(x = c(), y = c(), label = c())
for(i in 1:(length(boundaries)-1)) {
  point <- calc_map_xy_binned(vectorized_Sigmas, get_col_idx(vectorized_DUI, boundaries[i], boundaries[i+1]))
  plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, label = paste0("ABRP_",names(boundaries)[i+1])))
}

# repeat on Johnson et al. (2019)
vectorized_subject1 <- vectorized_Sigmas_johnson2019[1,]
ecdf <- sort(abs(vectorized_subject1))
boundaries <- quantile(ecdf, probs = seq(from = 0, to = 1, by = (1/4)))

for(i in 1:(length(boundaries)-1)) {
  point <- calc_map_xy_binned(vectorized_Sigmas_johnson2019, get_col_idx(vectorized_subject1, boundaries[i], boundaries[i+1]))
  plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, label = paste0("human_",names(boundaries)[i+1])))
}

p_quantiles <- ggplot(plot_df) +
  geom_point(aes(x = x, y = y, color = label), size = 3) +
  xlim(c(0, 1)) +
  ylim(c(-0.1, 1)) +
  xlab("avg. strength of associations (within hosts)") +
  ylab("avg. agreement between hosts")
print(p_quantiles)


























