library(ROL)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(abind)
library(matrixsampling)
library(stray)
library(driver)
library(lubridate)
library(stringr)

use_MAP <- TRUE

get_pairwise_combos <- function(P) {
  P^2 - P^2/2 - P/2
}

filter_taxa <- function(counts) {
  # filter to some minimum relative average abundance
  proportions <- as.matrix(counts)
  proportions <- apply(proportions, 2, function(x) x / sum(x))
  mean_rel_abundance <- rowMeans(proportions)
  mean_rel_abundance >= 0.001
}

## --------------------------------------------------------------------------------------------------------
##     Parse ABRP data
## --------------------------------------------------------------------------------------------------------

cat("Loading ABRP data...\n")

if(use_MAP) {
  Sigmas <- load_MAP_estimates(tax_level = "ASV", logratio = "clr")
  Sigmas <- lapply(Sigmas, function(x) cov2cor(x)) # JOHANNES, THANK YOU
} else {
  Sigmas <- load_full_posteriors(tax_level = "ASV", logratio = "clr")
  k <- dim(Sigmas[[1]])[3]
  for(i in 1:length(Sigmas)) {
    for(j in 1:k) {
      Sigmas[[i]][,,j] <- cov2cor(Sigmas[[i]][,,j])
    }
  }
}
n_hosts <- length(Sigmas)
P <- dim(Sigmas[[1]])[1]
pairwise_combos <- get_pairwise_combos(P)

# strip out the triangular portion of the correlation matrices, vectorize these, and put them in a matrix
# for later convenience
if(use_MAP) {
  vectorized_Sigmas <- matrix(NA, n_hosts, pairwise_combos)
  for(i in 1:n_hosts){
    vectorized_Sigmas[i,] <- Sigmas[[i]][upper.tri(Sigmas[[i]], diag = F)]
  }
} else {
  vectorized_Sigmas <- array(NA, dim = c(n_hosts, pairwise_combos, dim(Sigmas[[1]])[3]))
  for(i in 1:n_hosts){
    for(j in 1:k) {
      temp <- Sigmas[[i]][,,j]
      vectorized_Sigmas[i,,j] <- temp[upper.tri(temp, diag = F)]
    }
  }
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
##     Parse unfiltered Johnson et al. (2019) data; MAP estimates only
## --------------------------------------------------------------------------------------------------------

data_file <- file.path("input","fit_johnson.rds")
if(file.exists(data_file)) {
  vectorized_Sigmas_johnson2019 <- readRDS(data_file)
} else {
  # read count table
  counts <- read.table("input/johnson2019/taxonomy_counts_s.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # read sample ID to subject ID mapping
  mapping <- read.table("input/johnson2019/SampleID_map.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # filter to some minimum relative average abundance
  retain_idx <- filter_taxa(counts[,2:ncol(counts)])
  counts <- counts[retain_idx,]
  
  if(use_MAP) {
    depth <- 1
  } else {
    depth <- 20
  }
  
  D <- nrow(counts)
  n_interactions <- (D^2)/2 - D/2
  vectorized_Sigmas_johnson2019 <- array(NA, dim = c(length(unique(mapping$UserName)), n_interactions, depth))
  subjects <- unique(mapping$UserName)
  
  for(i in 1:length(subjects)) {
    subject <- subjects[i]
    cat("Evaluating subject:",subject,"\n")
    subject_samples <- mapping[mapping$UserName == subject,]$SampleID
    subject_days <- mapping[mapping$UserName == subject,]$StudyDayNo
    samples_to_pull <- colnames(counts)[which(colnames(counts) %in% subject_samples)]
    if(length(samples_to_pull) > 0) {
      subject_days <- which(subject_samples %in% samples_to_pull)
      subject_counts <- counts[,samples_to_pull] # omitting taxonomy
      # later: we probably want to this about removing taxa that are very rare within any individual
      
      # fit this with stray::basset
      Y <- as.matrix(subject_counts)
      X <- matrix(subject_days, 1, length(subject_days))
      
      alr_ys <- driver::alr((t(Y) + 0.5))
      alr_means <- colMeans(alr_ys)
      Theta <- function(X) matrix(alr_means, D-1, ncol(X))
      
      taxa_covariance <- get_Xi(D, total_variance = 1)
      
      rho <- calc_se_decay(min_correlation = 0.1, days_to_baseline = 1)
      Gamma <- function(X) {
        SE(X, sigma = 1, rho = rho, jitter = 1e-08)
      }
      
      if(use_MAP) {
        n_samples <- 0
        ret_mean <- TRUE
      } else {
        n_samples <- 20
        ret_mean <- FALSE
      }
        
      fit <- stray::basset(Y, X, taxa_covariance$upsilon, Theta, Gamma, taxa_covariance$Xi,
                           n_samples = n_samples, ret_mean = ret_mean,
                           b2 = 0.98, step_size = 0.004, eps_f = 1e-11, eps_g = 1e-05,
                           max_iter = 10000L, optim_method = "adam")
      fit.clr <- to_clr(fit)
      Sigmas <- fit.clr$Sigma
      for(k in 1:depth) {
        temp <- cov2cor(Sigmas[,,k])
        temp <- c(temp[upper.tri(temp, diag = F)])
        vectorized_Sigmas_johnson2019[i,,k] <- temp
      }
    } else {
      cat("\tSkipping subject",subject,"\n")
    }
  }

  # subject to complete cases (e.g. remove rows for the individuals with too few samples)
  include_rows <- which(!is.na(vectorized_Sigmas_johnson2019[,1,1]))
  vectorized_Sigmas_johnson2019 <- vectorized_Sigmas_johnson2019[include_rows,,,drop=F]
  saveRDS(vectorized_Sigmas_johnson2019, file = data_file)
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
##     Parse unfiltered Dethlefsen & Relman (2011) data; MAP estimates only
## --------------------------------------------------------------------------------------------------------

data_file <- file.path("input","fit_dethlefsen.rds")
if(file.exists(data_file)) {
  vectorized_Sigmas_dethlefsenrelman2011 <- readRDS(data_file)
} else {
  # read count table
  counts <- read.table("input/dethlefsen_relman2011/sd01.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  counts <- counts[,2:163]
  
  # filter low abundance taxa
  # counts <- counts[rowMeans(counts) > 20,]
  retain_idx <- filter_taxa(counts)
  counts <- counts[retain_idx,]
  
  # read the sampling schedule metadata
  metadata <- read.table("input/dethlefsen_relman2011/sampling_schedule.txt", header = FALSE, stringsAsFactors = FALSE)
  colnames(metadata) <- c("sample", "date")
  metadata$date <- as.Date(mdy(metadata$date))
  
  host_columns <- list(1:56, 57:108, 109:162) # individuals 1, 2, 3
  
  if(use_MAP) {
    depth <- 1
  } else {
    depth <- 20
  }
  
  D <- nrow(counts)
  n_interactions <- (D^2)/2 - D/2
  vectorized_Sigmas_dethlefsenrelman2011 <- array(NA, dim = c(length(host_columns), n_interactions, depth))
  
  for(subject in 1:length(host_columns)) {
    cat("Evaluating subject:",subject,"\n")
    subject_samples <- host_columns[[subject]]
    subject_dates <- metadata[metadata$sample %in% colnames(counts[subject_samples]),]$date
    # these appear to be pulled in order
    subject_days <- sapply(subject_dates, function(x) {
      difftime(x, subject_dates[1], units = "days")
    })
    if(length(subject_samples) > 0) {
      subject_counts <- counts[,subject_samples] # omitting taxonomy
      dim(subject_counts)
      # later: we probably want to this about removing taxa that are very rare within any individual
      
      # fit this with stray::basset
      Y <- as.matrix(subject_counts)
      X <- matrix(subject_days, 1, length(subject_days))
      
      alr_ys <- driver::alr((t(Y) + 0.5))
      alr_means <- colMeans(alr_ys)
      Theta <- function(X) matrix(alr_means, D-1, ncol(X))
      
      taxa_covariance <- get_Xi(D, total_variance = 1)
      
      rho <- calc_se_decay(min_correlation = 0.1, days_to_baseline = 7)
      Gamma <- function(X) {
        SE(X, sigma = 1, rho = rho, jitter = 1e-08)
      }
      
      if(use_MAP) {
        n_samples <- 0
        ret_mean <- TRUE
      } else {
        n_samples <- 20
        ret_mean <- FALSE
      }
      
      fit <- stray::basset(Y, X, taxa_covariance$upsilon, Theta, Gamma, taxa_covariance$Xi,
                           n_samples = n_samples, ret_mean = ret_mean,
                           b2 = 0.98, step_size = 0.004, eps_f = 1e-11, eps_g = 1e-05,
                           max_iter = 10000L, optim_method = "adam")
      
      fit.clr <- to_clr(fit)
      Sigmas <- fit.clr$Sigma
      for(k in 1:depth) {
        temp <- cov2cor(Sigmas[,,k])
        temp <- c(temp[upper.tri(temp, diag = F)])
        vectorized_Sigmas_dethlefsenrelman2011[subject,,k] <- temp
      }
    } else {
      cat("\tSkipping subject",subject,"\n")
    }
  }
  saveRDS(vectorized_Sigmas_dethlefsenrelman2011, file = data_file)
}

# # visualize the last fitted model for plausibility
# fit <- stray::basset(Y, X, taxa_covariance$upsilon, Theta, Gamma, taxa_covariance$Xi,
#                      n_samples = 500, ret_mean = FALSE)
# Eta <- predict(fit, subject_days, response = "Eta", iter = fit$iter)
# Eta <- alrInv_array(Eta, fit$D, 1)
# Eta <- clr_array(Eta, 1)
# lr_ys <- clr(t(fit$Y) + 0.5)
# 
# coord <- sample(1:D)[1]
# observations <- subject_days
# lr_tidy <- gather_array(lr_ys, "logratio_value", "timepoint", "logratio_coord")
# 
# no_samples <- dim(Eta)[3]
# posterior_samples <- gather_array(Eta[coord,,], "logratio_value", "observation", "sample_number")
# 
# # get quantiles
# post_quantiles <- posterior_samples %>%
#   group_by(observation) %>%
#   summarise(p2.5 = quantile(logratio_value, prob=0.025),
#            p5 = quantile(logratio_value, prob=0.05),
#            p10 = quantile(logratio_value, prob=0.1),
#            p25 = quantile(logratio_value, prob=0.25),
#            p50 = quantile(logratio_value, prob=0.5),
#            mean = mean(logratio_value),
#           p75 = quantile(logratio_value, prob=0.75),
#            p90 = quantile(logratio_value, prob=0.9),
#            p95 = quantile(logratio_value, prob=0.95),
#            p97.5 = quantile(logratio_value, prob=0.975)) %>%
#  ungroup()
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
##     Parse unfiltered Caporaso (2011) data
## --------------------------------------------------------------------------------------------------------

data_file <- file.path("input","fit_caporaso.rds")
if(file.exists(data_file)) {
  vectorized_Sigmas_caporaso2011 <- readRDS(data_file)
} else {
  # read count tables; OTUs all agree -- already checked
  counts_1 <- read.table("input/caporaso2011/F4_feces_L6.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  counts_1 <- counts_1[,2:ncol(counts_1)]
  
  counts_2 <- read.table("input/caporaso2011/M3_feces_L6.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  counts_2 <- counts_2[,2:ncol(counts_2)]
  
  # first row are the sample day identifiers
  counts <- cbind(counts_1, counts_2)
  
  # filter low abundance taxa
  #counts <- counts[rowMeans(counts) > 20,]
  retain_idx <- filter_taxa(counts[2:nrow(counts),])
  counts <- counts[c(TRUE, retain_idx),]
  
  host_columns <- list(1:ncol(counts_1), (ncol(counts_1)+1):ncol(counts)) # individuals 1, 2
  pairwise_combos <- get_pairwise_combos(nrow(counts) - 1)
  
  if(use_MAP) {
    depth <- 1
  } else {
    depth <- 20
  }
  
  D <- nrow(counts) - 1
  n_interactions <- (D^2)/2 - D/2
  vectorized_Sigmas_caporaso2011 <- array(NA, dim = c(length(host_columns), n_interactions, depth))
  
  for(subject in 1:length(host_columns)) {
    cat("Evaluating subject:",subject,"\n")
    subject_samples <- host_columns[[subject]]
    subject_dates <- counts[1,host_columns[[subject]]]
    if(length(subject_samples) > 0) {
      subject_counts <- counts[2:nrow(counts),subject_samples] # omitting taxonomy
      # later: we probably want to this about removing taxa that are very rare within any individual
      
      # fit this with stray::basset
      Y <- as.matrix(subject_counts)
      X <- matrix(subject_dates, 1, length(subject_dates))
      
      alr_ys <- driver::alr((t(Y) + 0.5))
      alr_means <- colMeans(alr_ys)
      Theta <- function(X) matrix(alr_means, D-1, ncol(X))
      
      taxa_covariance <- get_Xi(D, total_variance = 1)
      
      rho <- calc_se_decay(min_correlation = 0.1, days_to_baseline = 7)
      Gamma <- function(X) {
        SE(X, sigma = 1, rho = rho, jitter = 1e-08)
      }
      
      if(use_MAP) {
        n_samples <- 0
        ret_mean <- TRUE
      } else {
        n_samples <- 20
        ret_mean <- FALSE
      }
      
      # full data set
      fit <- stray::basset(Y, X, taxa_covariance$upsilon, Theta, Gamma, taxa_covariance$Xi,
                           n_samples = n_samples, ret_mean = ret_mean,
                           b2 = 0.98, step_size = 0.004, eps_f = 1e-11, eps_g = 1e-05,
                           max_iter = 10000L, optim_method = "adam")
      
      fit.clr <- to_clr(fit)
      Sigmas <- fit.clr$Sigma
      for(k in 1:depth) {
        temp <- cov2cor(Sigmas[,,k])
        temp <- c(temp[upper.tri(temp, diag = F)])
        vectorized_Sigmas_caporaso2011[subject,,k] <- temp
      }
    } else {
      cat("\tSkipping subject",subject,"\n")
    }
  }
  saveRDS(vectorized_Sigmas_caporaso2011, file = data_file)
}

# n_subsets <- 20
# vectorized_Sigmas_caporaso2011_subsetted <- array(NA, dim = c(length(host_columns), n_interactions, depth, n_subsets))
# 
# for(m in 1:n_subsets) {
#   # subset the Caporaso data to the same number of samples as the Johnson et al. data
#   host_columns <- list(1:ncol(counts_1), (ncol(counts_1)+1):ncol(counts)) # individuals 1, 2
#   # subset CONSEQUTIVE in time
#   # offset1 <- floor(runif(1, min = 1, max = length(host_columns[[1]]) - 17))
#   # offset2 <- floor(runif(1, min = 1, max = length(host_columns[[2]]) - 17))
#   # host_columns[[1]] <- host_columns[[1]][offset1:(offset1+16)]
#   # host_columns[[2]] <- host_columns[[2]][offset2:(offset2+16)]
#   # subset GAPPED in time
#   host_columns[[1]] <- sort(sample(host_columns[[1]])[1:17])
#   host_columns[[2]] <- sort(sample(host_columns[[2]])[1:17])
#   
#   for(subject in 1:length(host_columns)) {
#     cat("Evaluating subject:",subject,"\n")
#     subject_samples <- host_columns[[subject]]
#     subject_dates <- counts[1,host_columns[[subject]]]
#     if(length(subject_samples) > 0) {
#       subject_counts <- counts[2:nrow(counts),subject_samples] # omitting taxonomy
#       # later: we probably want to this about removing taxa that are very rare within any individual
#       
#       # fit this with stray::basset
#       Y <- as.matrix(subject_counts)
#       X <- matrix(subject_dates, 1, length(subject_dates))
#       
#       alr_ys <- driver::alr((t(Y) + 0.5))
#       alr_means <- colMeans(alr_ys)
#       Theta <- function(X) matrix(alr_means, D-1, ncol(X))
#       
#       taxa_covariance <- get_Xi(D, total_variance = 1)
#       
#       rho <- calc_se_decay(min_correlation = 0.1, days_to_baseline = 7)
#       Gamma <- function(X) {
#         SE(X, sigma = 1, rho = rho, jitter = 1e-08)
#       }
#       
#       if(use_MAP) {
#         n_samples <- 0
#         ret_mean <- TRUE
#       } else {
#         n_samples <- 20
#         ret_mean <- FALSE
#       }
#       
#       # full data set
#       fit <- stray::basset(Y, X, taxa_covariance$upsilon, Theta, Gamma, taxa_covariance$Xi,
#                            n_samples = n_samples, ret_mean = ret_mean,
#                            b2 = 0.98, step_size = 0.004, eps_f = 1e-11, eps_g = 1e-05,
#                            max_iter = 10000L, optim_method = "adam")
#       
#       fit.clr <- to_clr(fit)
#       Sigmas <- fit.clr$Sigma
#       for(k in 1:depth) {
#         temp <- cov2cor(Sigmas[,,k])
#         temp <- c(temp[upper.tri(temp, diag = F)])
#         vectorized_Sigmas_caporaso2011_subsetted[subject,,k,m] <- temp
#       }
#     } else {
#       cat("\tSkipping subject",subject,"\n")
#     }
#   }
# }

# # visualize the last fitted model for plausibility
# fit <- stray::basset(Y, X, taxa_covariance$upsilon, Theta, Gamma, taxa_covariance$Xi,
#                      n_samples = 500, ret_mean = FALSE)
# #                     max_iter = 10000L, optim_method = "adam")
# Eta <- predict(fit, subject_dates, response = "Eta", iter = fit$iter)
# Eta <- alrInv_array(Eta, fit$D, 1)
# Eta <- clr_array(Eta, 1)
# lr_ys <- clr(t(fit$Y) + 0.5)
# 
# coord <- sample(1:D)[1]
# observations <- subject_dates
# lr_tidy <- gather_array(lr_ys, "logratio_value", "timepoint", "logratio_coord")
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
##     Parse David (2014) data
## --------------------------------------------------------------------------------------------------------

data_file <- file.path("input","fit_david.rds")
if(file.exists(data_file)) {
  vectorized_Sigmas_david2014 <- readRDS(data_file)
} else {
  counts <- read.table("input/david2014/otu.table.ggref", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  metadata <- read.table("input/david2014/13059_2013_3286_MOESM18_ESM.csv", sep = ",", header = TRUE)
  
  sample_labels <- metadata$X
  subj_labels <- metadata$AGE # subject A is 26, subject B is 36
  day_labels <- metadata$COLLECTION_DAY
  subj_labels[subj_labels == 26] <- "A"
  subj_labels[subj_labels == 36] <- "B"
  # strip character from the count column names
  sample_labels <- sapply(sample_labels, function(x) {
    str_replace(x, "\\.\\d+", "")
  })
  sample_labels <- unname(sample_labels)

  # remove saliva samples
  keep_idx <- which(!sapply(sample_labels, function(x) {
    str_detect(x, "Saliva")
  }))
  subj_labels <- subj_labels[keep_idx]
  sample_labels <- sample_labels[keep_idx]
  day_labels <- day_labels[keep_idx]
  
  # include columns in counts that are in sample_labels
  counts <- counts[,colnames(counts) %in% sample_labels]
  keep_idx <- sample_labels %in% colnames(counts)
  sample_labels <- sample_labels[keep_idx]
  subj_labels <- subj_labels[keep_idx]
  day_labels <- day_labels[keep_idx]

  subj_A_idx <- c()
  subj_A_days <- c()
  subj_B_idx <- c()
  subj_B_days <- c()
  for(i in 1:length(colnames(counts))) {
    idx <- which(sample_labels == colnames(counts)[i])
    if(subj_labels[idx] == "A") {
      subj_A_idx <- c(subj_A_idx, i)
      subj_A_days <- c(subj_A_days, day_labels[idx])
    } else {
      subj_B_idx <- c(subj_B_idx, i)
      subj_B_days <- c(subj_B_days, day_labels[idx])
    }
  }
  
  subj_A_counts <- counts[,subj_A_idx]
  subj_B_counts <- counts[,subj_B_idx]

  reorder <- order(subj_A_days)
  subj_A_days <- subj_A_days[reorder]
  subj_A_counts <- subj_A_counts[,reorder]

  reorder <- order(subj_B_days)
  subj_B_days <- subj_B_days[reorder]
  subj_B_counts <- subj_B_counts[,reorder]
  
  counts <- cbind(subj_A_counts, subj_B_counts)
  # filter low abundance taxa
  # counts <- counts[rowMeans(counts) > 100,]
  retain_idx <- filter_taxa(counts)
  counts <- counts[retain_idx,]
  
  host_columns <- list(1:length(subj_A_days), (length(subj_A_days) + 1):ncol(counts)) # individuals 1, 2
  host_dates <- list(subj_A_days, subj_B_days)

  # for testing
  #host_columns[[1]] <- host_columns[[1]][1:50]
  #host_columns[[2]] <- host_columns[[2]][1:50]
  #host_dates[[1]] <- host_dates[[1]][1:50]
  #host_dates[[2]] <- host_dates[[2]][1:50]
    
  if(use_MAP) {
    depth <- 1
  } else {
    depth <- 20
  }
  
  D <- nrow(counts)
  n_interactions <- (D^2)/2 - D/2
  vectorized_Sigmas_david2014 <- array(NA, dim = c(length(host_columns), n_interactions, depth))
  
  for(subject in 1:length(host_columns)) {
    cat("Evaluating subject:",subject,"\n")
    subject_samples <- host_columns[[subject]]
    subject_dates <- host_dates[[subject]] + 1
    if(length(subject_samples) > 0) {
      subject_counts <- counts[,subject_samples] # omitting taxonomy
      # later: we probably want to this about removing taxa that are very rare within any individual
      
      # fit this with stray::basset
      Y <- as.matrix(subject_counts)
      X <- matrix(subject_dates, 1, length(subject_dates))
      
      alr_ys <- driver::alr((t(Y) + 0.5))
      alr_means <- colMeans(alr_ys)
      Theta <- function(X) matrix(alr_means, D-1, ncol(X))
      
      taxa_covariance <- get_Xi(D, total_variance = 1)
      
      rho <- calc_se_decay(min_correlation = 0.1, days_to_baseline = 7)
      Gamma <- function(X) {
        SE(X, sigma = 1, rho = rho, jitter = 1e-08)
      }
      
      if(use_MAP) {
        n_samples <- 0
        ret_mean <- TRUE
      } else {
        n_samples <- 20
        ret_mean <- FALSE
      }
      
      # full data set
      fit <- stray::basset(Y, X, taxa_covariance$upsilon, Theta, Gamma, taxa_covariance$Xi,
                           n_samples = n_samples, ret_mean = ret_mean,
                           b2 = 0.98, step_size = 0.004, eps_f = 1e-11, eps_g = 1e-05,
                           max_iter = 10000L, optim_method = "adam")
      
      fit.clr <- to_clr(fit)
      Sigmas <- fit.clr$Sigma
      for(k in 1:depth) {
        temp <- cov2cor(Sigmas[,,k])
        temp <- c(temp[upper.tri(temp, diag = F)])
        vectorized_Sigmas_david2014[subject,,k] <- temp
      }
    } else {
      cat("\tSkipping subject",subject,"\n")
    }
  }
  saveRDS(vectorized_Sigmas_david2014, file = data_file)
}

# # visualize the last fitted model for plausibility
# fit <- stray::basset(Y, X, taxa_covariance$upsilon, Theta, Gamma, taxa_covariance$Xi,
#                     n_samples = 500, ret_mean = FALSE)
# #                     max_iter = 10000L, optim_method = "adam")
# Eta <- predict(fit, subject_dates, response = "Eta", iter = fit$iter)
# Eta <- alrInv_array(Eta, fit$D, 1)
# Eta <- clr_array(Eta, 1)
# lr_ys <- clr(t(fit$Y) + 0.5)
# 
# coord <- sample(1:D)[1]
# observations <- subject_dates
# lr_tidy <- gather_array(lr_ys, "logratio_value", "timepoint", "logratio_coord")
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
##     Visualize ABRP social group vs. human data outgroups
## --------------------------------------------------------------------------------------------------------

calc_map_xy <- function(vectorized_Sigmas) {
  within_score <- mean(apply(abs(vectorized_Sigmas), 1, mean))
  between_score <- mean(apply(combn(1:nrow(vectorized_Sigmas), m = 2), 2, function(x) {
    cor(vectorized_Sigmas[x[1],], vectorized_Sigmas[x[2],])
  }))
  return(list(x = within_score, y = between_score))
}

if(use_MAP) {
  plot_df <- data.frame(x = c(), y = c(), group = c())
  point <- calc_map_xy(vectorized_Sigmas_johnson2019[,,1])
  plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = "Johnson et al. (2019)"))
  point <- calc_map_xy(vectorized_Sigmas_dethlefsenrelman2011[,,1])
  plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = "Dethlefsen & Relman (2011)"))
  point <- calc_map_xy(vectorized_Sigmas_caporaso2011[,,1])
  plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = "Caporaso et al. (2011)"))
  point <- calc_map_xy(vectorized_Sigmas_david2014[,,1])
  plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = "David et al. (2014)"))
  # for(m in 1:n_subsets) {
  #   point <- calc_map_xy(vectorized_Sigmas_caporaso2011_subsetted[,,1,m])
  #   plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = "Caporaso et al. (2011), subsetted"))
  # }
  for(group in unlist(unique(labels))) {
    group_Sigmas <- vectorized_Sigmas[unname(which(labels == group)),]
    point <- calc_map_xy(group_Sigmas)
    plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = paste0("ABRP group ",group)))
  }
} else {
  plot_df <- data.frame(x = c(), y = c(), group = c())
  # Johnson et al. (2019)
  for(j in 1:dim(vectorized_Sigmas_johnson2019)[3]) {
    point <- calc_map_xy(vectorized_Sigmas_johnson2019[,,j])
    plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = "Johnson et al. (2019)"))
  }
  # Dethlefsen & Relman (2011)
  for(j in 1:dim(vectorized_Sigmas_dethlefsenrelman2011)[3]) {
    point <- calc_map_xy(vectorized_Sigmas_dethlefsenrelman2011[,,j])
    plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = "Dethlefsen, Relman (2011)"))
  }
  # Caporaso et al. (2011)
  for(j in 1:dim(vectorized_Sigmas_caporaso2011)[3]) {
    point <- calc_map_xy(vectorized_Sigmas_caporaso2011[,,j])
    plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = "Caporaso et al. (2011)"))
  }
  # David et al. (2014)
  for(j in 1:dim(vectorized_Sigmas_david2014)[3]) {
    point <- calc_map_xy(vectorized_Sigmas_david2014[,,j])
    plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = "David et al. (2014)"))
  }
  # ABRP
  for(group in unlist(unique(labels))) {
    group_Sigmas <- vectorized_Sigmas[unname(which(labels == group)),,]
    for(j in 1:k) {
      point <- calc_map_xy(group_Sigmas[,,j])
      plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = paste0("ABRP group ",group)))
    }
  }
}

p_combo <- ggplot(plot_df) +
  geom_point(aes(x = y, y = x, color = group), size = 4, stroke = 1) +
  scale_color_manual(values = c(
                                "#3486eb", # JOHNSON
                                "#32d1bf", # DETHLEFSEN
                                "#87d132", # CAPORASO
                                "#bd34eb", # DAVID
                                "#eb7434", # ABRP
                                "#eb8034", # ABRP
                                "#eb8f34", # ABRP
                                "#eba234", # ABRP
                                "#ebb434"  # ABRP
                                )) +
  xlim(c(0,1)) +
  ylim(c(-0.1,1)) +
  xlab("avg. agreement between hosts") +
  ylab("avg. strength of associations (within hosts)")
# print(p_combo)
ggsave("map_human_datasets.png", p_combo, units = "in", dpi = 100, height = 8, width = 10)

## --------------------------------------------------------------------------------------------------------
##     Visualize on "map" as quantiles
## --------------------------------------------------------------------------------------------------------

# vectorized_DUI <- vectorized_Sigmas[which(names(Sigmas) == "DUI"),]
# ecdf <- sort(abs(vectorized_DUI))
# # separate into quantiles
# boundaries <- quantile(ecdf, probs = seq(from = 0, to = 1, by = 1))
# 
# # bin a set of pairwise relationships between microbes
# get_col_idx <- function(vec, lower_boundary, upper_boundary) {
#   which(abs(vec) >= lower_boundary & abs(vec) < upper_boundary)
# }
# 
# # col_idx are the column indices of pairs of microbes to consider here
# calc_map_xy_binned <- function(vectorized_Sigmas, col_idx) {
#   within_score <- mean(apply(abs(vectorized_Sigmas[,col_idx]), 1, mean))
#   between_score <- mean(apply(combn(1:nrow(vectorized_Sigmas), m = 2), 2, function(x) {
#     cor(vectorized_Sigmas[x[1],col_idx], vectorized_Sigmas[x[2],col_idx])
#   }))
#   return(list(x = within_score, y = between_score))
# }
# 
# # add in baboon quantiles
# plot_df <- data.frame(x = c(), y = c(), label = c())
# for(i in 1:(length(boundaries)-1)) {
#   point <- calc_map_xy_binned(vectorized_Sigmas, get_col_idx(vectorized_DUI, boundaries[i], boundaries[i+1]))
#   plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, label = paste0("ABRP_",names(boundaries)[i+1])))
# }
# 
# # repeat on Johnson et al. (2019)
# # vectorized_subject1 <- vectorized_Sigmas_johnson2019[1,]
# # ecdf <- sort(abs(vectorized_subject1))
# # boundaries <- quantile(ecdf, probs = seq(from = 0, to = 1, by = (1/4)))
# # for(i in 1:(length(boundaries)-1)) {
# #   point <- calc_map_xy_binned(vectorized_Sigmas_johnson2019, get_col_idx(vectorized_subject1, boundaries[i], boundaries[i+1]))
# #   plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, label = paste0("human_",names(boundaries)[i+1])))
# # }
# 
# p_quantiles <- ggplot(plot_df) +
#   geom_point(aes(x = y, y = x, color = label), size = 4) +
#   xlim(c(0, 1)) +
#   ylim(c(-0.1, 1)) +
#   xlab("avg. agreement between hosts") +
#   ylab("avg. strength of associations (within hosts)")
# print(p_quantiles)

