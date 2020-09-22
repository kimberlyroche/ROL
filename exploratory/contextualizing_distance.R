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

# counts is taxa x samples
# host_columns is a list (length = num. hosts) of host column indices in counts
# host_dates is a list (length = num. hosts) of host sample dates associated with the columns in counts
# depth is the number of posterior samples to draw (depth = 1 is MAP estimation)
fit_model <- function(counts, host_columns, host_dates, depth = 1) {
  D <- nrow(counts)
  n_interactions <- (D^2)/2 - D/2
  vectorized_Sigmas <- array(NA, dim = c(length(host_columns), n_interactions, depth))
  
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
      
      if(depth == 1) {
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
        vectorized_Sigmas[subject,,k] <- temp
      }
    } else {
      cat("\tSkipping subject",subject,"\n")
    }
  }
  return(vectorized_Sigmas)
}

# TO DO: general purpose fit visualization function
# visualize the last fitted model for plausibility
# fit <- stray::basset(Y, X, taxa_covariance$upsilon, Theta, Gamma, taxa_covariance$Xi,
#                     n_samples = 500, ret_mean = FALSE)
# #                     max_iter = 10000L, optim_method = "adam")
# Eta <- predict(fit, min(subject_dates):max(subject_dates), response = "Eta", iter = fit$iter)
# Eta <- alrInv_array(Eta, fit$D, 1)
# Eta <- clr_array(Eta, 1)
# lr_ys <- clr(t(fit$Y) + 0.5)
# 
# coord <- sample(1:D)[1]
# lr_coord <- data.frame(x = subject_dates, y = lr_ys[,coord])
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
#   geom_point(data = lr_coord, aes(x = x, y = y), alpha=0.5) +
#   theme_minimal() +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(angle=45)) +
#   ylab("LR coord")
# show(p)

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

  subjects <- unique(mapping$UserName)  
  host_columns <- list()
  host_dates <- list()
  for(i in 1:length(subjects)) {
    subject <- subjects[i]
    subject_sample_IDs <- mapping[mapping$UserName == subject,]$SampleID
    subject_days <- mapping[mapping$UserName == subject,]$StudyDayNo
    subject_columns <- which(colnames(counts) %in% subject_sample_IDs)
    subject_days <- which(subject_sample_IDs %in% colnames(counts)[subject_columns])
    if(length(subject_days) > 0) {
      subject_days <- subject_days - min(subject_days)
      idx <- length(host_columns) + 1
      host_columns[[idx]] <- subject_columns
      host_dates[[idx]] <- subject_days
    }
  }
  
  vectorized_Sigmas_johnson2019 <- fit_model(counts, host_columns, host_dates)
  saveRDS(vectorized_Sigmas_johnson2019, file = data_file)
}

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
  host_dates <- list()
  for(i in 1:3) {
    subject_dates <- metadata[metadata$sample %in% colnames(counts[host_columns[[i]]]),]$date
    host_dates[[i]] <- sapply(subject_dates, function(x) {
      difftime(x, subject_dates[1], units = "days")
    })
  }
  
  vectorized_Sigmas_dethlefsenrelman2011 <- fit_model(counts, host_columns, host_dates)
  saveRDS(vectorized_Sigmas_dethlefsenrelman2011, file = data_file)
}

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
  host_dates <- list()
  for(i in 1:2) {
    host_dates[[i]] <- unname(unlist(counts[1,host_columns[[i]]]))
  }
  
  vectorized_Sigmas_caporaso2011 <- fit_model(counts, host_columns, host_dates)
  saveRDS(vectorized_Sigmas_caporaso2011, file = data_file)
}

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

  vectorized_Sigmas_david2014 <- fit_model(counts, host_columns, host_dates)
  saveRDS(vectorized_Sigmas_david2014, file = data_file)
}

## --------------------------------------------------------------------------------------------------------
##     Parse DIABIMMUNE (infants)
## --------------------------------------------------------------------------------------------------------

data_file <- file.path("input","fit_DIABIMMUNE.rds")
if(file.exists(data_file)) {
  vectorized_Sigmas_DIABIMMUNE <- readRDS(data_file)
} else {
  load("input/DIABIMMUNE/diabimmune_karelia_16s_data.rdata")
  load("input/DIABIMMUNE/DIABIMMUNE_Karelia_metadata.RData")
  
  counts <- t(data_16s)
  rm(data_16s)
  
  # these aren't counts but relative abundances
  # let's dream and scale them into counts
  for(i in 1:ncol(counts)) {
    counts[,i] <- rmultinom(1, size = rpois(1, 10000), prob = counts[,i])
  }
  retain_idx <- filter_taxa(counts)
  counts <- counts[retain_idx,]
  
  # use kids with at least 15 samples
  use_subjects <- names(which(table(metadata$subjectID) >= 15))
  #metadata <- metadata[metadata$subjectID %in% use_subjects,]

  host_columns <- list()
  host_dates <- list()
  for(i in 1:ncol(counts)) {
    sample_ID <- colnames(counts)[i]
    metadata_idx <- which(metadata$SampleID == sample_ID)
    subject_ID <- metadata$subjectID[metadata_idx]
    if(subject_ID %in% use_subjects) {
      age_at_collection <- metadata$age_at_collection[metadata_idx]
      if(subject_ID %in% names(host_columns)) {
        host_columns[[subject_ID]] <- c(host_columns[[subject_ID]], i)
        host_dates[[subject_ID]] <- c(host_dates[[subject_ID]], age_at_collection)
      } else {
        host_columns[[subject_ID]] <- c(i)
        host_dates[[subject_ID]] <- age_at_collection
      }
    }
  }
  
  for(i in 1:length(host_dates)) {
    baseline_day <- min(host_dates[[i]])
    sample_days <- sapply(host_dates[[i]], function(x) x - baseline_day) + 1
    host_dates[[i]] <- sample_days
  }
  
  vectorized_Sigmas_DIABIMMUNE <- fit_model(counts, host_columns, host_dates)
  saveRDS(vectorized_Sigmas_DIABIMMUNE, file = data_file)
}

## --------------------------------------------------------------------------------------------------------
##     Parse Grossart lakes data (Germany)
## --------------------------------------------------------------------------------------------------------

data_file <- file.path("input","fit_grossart.rds")
if(file.exists(data_file)) {
  vectorized_Sigmas_grossart <- readRDS(data_file)
} else {
  counts <- readRDS("input/grossart_lakes/data.rds")
  counts <- otu_table(counts)@.Data
  retain_idx <- filter_taxa(counts)
  counts <- counts[retain_idx,]
  
  metadata <- read.table("input/grossart_lakes/945_20190102-074005.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  # look at the sample numbers per unique combo of depth and filter poresize for each lake
  # we'll hand pick similar conditions with large sample numbers!
  # temp <- metadata %>%
  #   group_by(description, depth, filter_poresize) %>%
  #   tally()
  # print(as_tibble(temp), n = 100)
  
  conditions <- list(c("freshwater metagenome, Lake Breiter Luzin", "0-5", "0.2micron"),
                     c("freshwater metagenome, Lake Grosse Fuchskuhle", "0-2", "0.2micron"),
                     c("freshwater metagenome, Lake Melzer", "0-1", "0.2micron"),
                     c("freshwater metagenome, Lake Stechlin", "0-10", "0.2micron"),
                     c("freshwater metagenome, Lake Tiefwaren", "0-10", "0.2micron"))
  
  # cycle through and collect the samples for all of these
  host_columns <- list()
  host_dates <- list()
  for(i in 1:length(conditions)) {
    focal_md <- metadata[metadata$description == conditions[[i]][1] &
                           metadata$depth == conditions[[i]][2] &
                           metadata$filter_poresize == conditions[[i]][3],]
    sample_names <- focal_md$sample_name
    timestamps <- focal_md$collection_timestamp
    reorder <- order(timestamps)
    sample_names <- sample_names[reorder]
    timestamps <- timestamps[reorder]
    
    # it looks like some sample names aren't in the count table?
    matched_sample_names <- sample_names %in% colnames(counts)
    sample_names <- sample_names[matched_sample_names]
    timestamps <- timestamps[matched_sample_names]
    
    baseline_date <- timestamps[1]
    days_vector <- unname(sapply(timestamps, function(x) difftime(as.Date(x), as.Date(baseline_date), units = "days")) + 1)
    host_columns[[i]] <- which(colnames(counts) %in% sample_names)
    host_dates[[i]] <- days_vector
  }
  
  vectorized_Sigmas_grossart <- fit_model(counts, host_columns, host_dates)
  saveRDS(vectorized_Sigmas_grossart, file = data_file)
}

## --------------------------------------------------------------------------------------------------------
##     Parse McMahon lakes data (Wisconsin, USA)
## --------------------------------------------------------------------------------------------------------

data_file_1 <- file.path("input","fit_mcmahon_1.rds")
data_file_2 <- file.path("input","fit_mcmahon_1.rds")
if(file.exists(data_file_1) & file.exists(data_file_2)) {
  vectorized_Sigmas_mcmahon_E <- readRDS(data_file_1)
  vectorized_Sigmas_mcmahon_H <- readRDS(data_file_2)
} else {
  counts <- readRDS("input/mcmahon_lakes/data.rds")
  counts <- otu_table(counts)@.Data
  retain_idx <- filter_taxa(counts)
  counts <- counts[retain_idx,]
  
  metadata <- read.table("input/mcmahon_lakes/1288_20180418-110149.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

  # lake short name to list index lookup
  # lake_names <- c("MA", "TB", "CB", "NS", "WSB", "SSB", "HK", "NSB")
  metadata$lake <- NA
  metadata$layer <- NA
  # append more useful columns to the metadata
  for(i in 1:nrow(metadata)) {
    matches <- str_match_all(metadata$description[i], regex("^freshwater metagenome (\\D+)(\\d+)(\\D+)(\\d+)"))
    lake_name <- matches[[1]][1,2]
    if(lake_name != "NSb") {
      layer <- str_match_all(lake_name, regex("[H|E]$"))
      if(nrow(layer[[1]]) > 0) {
        if(layer[[1]][1,1] == "H") {
          # hypolimnion
          metadata$lake[i] <- substr(lake_name, 1, str_length(lake_name)-1)
          metadata$layer[i] <- "H"
        } else {
          # epilimnion
          metadata$lake[i] <- substr(lake_name, 1, str_length(lake_name)-1)
          metadata$layer[i] <- "E"
        }
      } else {
        metadata$lake[i] <- lake_name
      }
    }
  }
  
  host_columns_E <- list()
  host_columns_H <- list()
  host_dates_E <- list()
  host_dates_H <- list()
  for(i in 1:ncol(counts)) {
    sample_id <- colnames(counts)[i]
    metadata_idx <- which(metadata$sample_name == sample_id)
    lake_name <- metadata$lake[metadata_idx]
    layer <- metadata$layer[metadata_idx]
    if(!is.na(layer) & !is.na(lake_name)) { # "NSb", the singleton
      sample_date <- as.Date(metadata$collection_timestamp[metadata_idx])
      if(layer == "E") {
        if(lake_name %in% names(host_columns_E)) {
          host_columns_E[[lake_name]] <- c(host_columns_E[[lake_name]], i)
          host_dates_E[[lake_name]] <- c(host_dates_E[[lake_name]], sample_date)
        } else {
          host_columns_E[[lake_name]] <- c(i)
          host_dates_E[[lake_name]] <- sample_date
        }
      } else {
        if(lake_name %in% names(host_columns_H)) {
          host_columns_H[[lake_name]] <- c(host_columns_H[[lake_name]], i)
          host_dates_H[[lake_name]] <- c(host_dates_H[[lake_name]], sample_date)
        } else {
          host_columns_H[[lake_name]] <- c(i)
          host_dates_H[[lake_name]] <- sample_date
        }
      }
    }
  }

  # convert dates to days
  for(i in 1:length(host_dates_E)) {
    baseline_date <- min(host_dates_E[[i]])
    sample_days <- sapply(host_dates_E[[i]], function(x) difftime(x, baseline_date, units = "days")) + 1
    host_dates_E[[i]] <- sample_days
  }
  for(i in 1:length(host_dates_H)) {
    baseline_date <- min(host_dates_H[[i]])
    sample_days <- sapply(host_dates_H[[i]], function(x) difftime(x, baseline_date, units = "days")) + 1
    host_dates_H[[i]] <- sample_days
  }
  
  vectorized_Sigmas_mcmahon_E <- fit_model(counts, host_columns_E, host_dates_E)
  saveRDS(vectorized_Sigmas_mcmahon_E, file = data_file_1)
  vectorized_Sigmas_mcmahon_H <- fit_model(counts, host_columns_H, host_dates_H)
  saveRDS(vectorized_Sigmas_mcmahon_H, file = data_file_2)
}

## --------------------------------------------------------------------------------------------------------
##     Visualize ABRP social group vs. human data outgroups
## --------------------------------------------------------------------------------------------------------

# note: x and y are swapped on the rendered plot!
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
  point <- calc_map_xy(vectorized_Sigmas_DIABIMMUNE[,,1])
  plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = "DIABIMMUNE"))
  point <- calc_map_xy(vectorized_Sigmas_grossart[,,1])
  plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = "Grossart lakes data"))
  point <- calc_map_xy(vectorized_Sigmas_mcmahon_E[,,1])
  plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = "McMahon lakes data (epilimnion)"))
  point <- calc_map_xy(vectorized_Sigmas_mcmahon_H[,,1])
  plot_df <- rbind(plot_df, data.frame(x = point$x, y = point$y, group = "McMahon lakes data (hypolimnion)"))
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
                                "#ebb434", # ABRP
                                "#87d132", # CAPORASO
                                "#32d1bf", # DAVID
                                "#3486eb", # DETHLEFSEN
                                "#ffa7b6", # DIABIMMUNE
                                "#999999", # GROSSART
                                "#bd34eb", # JOHNSON
                                "#aaaaaa", # MCMAHON (epilimnion; shallow water)
                                "#bbbbbb"  # MCMAHON (hypolimnion; deep water)
  )) +
  xlim(c(0,1)) +
  ylim(c(-0.1,1)) +
  xlab("avg. agreement between hosts") +
  ylab("avg. strength of associations (within hosts)")
print(p_combo)
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

