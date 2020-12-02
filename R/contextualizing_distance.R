#' Summarizes a matrix of correlations (rows are hosts, columns are taxon-taxon pairs) in 2D
#' 
#' @param heatmap vector or matrix (rows as hosts) of taxon-taxon correlations across hosts
#' @return NULL
#' @export
summarize_all_pairs_2D <- function(vectorized_Sigmas) {
  # calculate average correlation between vectors for all pairs of hosts
  x_coord <- mean(apply(combn(1:nrow(vectorized_Sigmas), m = 2), 2, function(x) {
    cor(vectorized_Sigmas[x[1],], vectorized_Sigmas[x[2],])
  }))
  # calculate mean absolute correlation strength, typically low
  y_coord <- mean(apply(abs(vectorized_Sigmas), 2, mean))
  return(list(x = x_coord, y = y_coord))
}

#' Utility function to calculate the number of unique combinations of P features
#' 
#' @param P feature number
#' @return number of unique combinations of features
#' @export
n_pairwise_combos <- function(P) {
  P^2/2 - P/2
}

#' Given a count table, label taxa with at least 0.1% relative abundance
#' 
#' @param counts count table (taxa as rows, samples as columns)
#' @return number of unique combinations of features
#' @export
filter_taxa_2D_summary <- function(counts) {
  # filter to some minimum relative average abundance
  proportions <- as.matrix(counts)
  proportions <- apply(proportions, 2, function(x) x / sum(x))
  mean_rel_abundance <- rowMeans(proportions)
  mean_rel_abundance >= 0.001
}

#' Parse correlation matrix for Amboseli data
#' 
#' @param heatmap vector or matrix (rows as hosts) of taxon-taxon correlations across hosts
#' @return NULL
#' @details
#' Note: This version is using CLR correlation (not proportionality).
#' @import phyloseq
#' @export
load_ABRP_2D_summary <- function(MAP = TRUE) {
  if(MAP) {
    Sigmas <- load_MAP_estimates(tax_level = "ASV", logratio = "clr")
    Sigmas <- lapply(Sigmas, function(x) cov2cor(x))
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
  P <- dim(Sigmas[[1]])[1] # number of taxa
  n_pairwise_combos <- n_pairwise_combos(P)
  
  # Strip out the triangular portion of the correlation matrices, vectorize these, and put them
  # in a matrix (hosts x taxon-taxon pairs)
  if(MAP) {
    vectorized_Sigmas <- matrix(NA, n_hosts, n_pairwise_combos)
    for(i in 1:n_hosts){
      vectorized_Sigmas[i,] <- Sigmas[[i]][upper.tri(Sigmas[[i]], diag = FALSE)]
    }
  } else {
    vectorized_Sigmas <- array(NA, dim = c(n_hosts, n_pairwise_combos, dim(Sigmas[[1]])[3]))
    for(i in 1:n_hosts){
      for(j in 1:k) {
        temp <- Sigmas[[i]][,,j]
        vectorized_Sigmas[i,,j] <- temp[upper.tri(temp, diag = FALSE)]
      }
    }
  }
  
  # Get host primary social group assignments. Pull these from metadata of
  # phyloseq object.
  data <- load_data(tax_level = "ASV")
  labels <- list()
  for(host in names(Sigmas)) {
    host <<- host
    subdata <- subset_samples(data, sname == host)
    metadata <- sample_data(subdata)
    labels[[host]] <- names(which(table(metadata[["grp"]]) == max(table(metadata[["grp"]]))))[1]
  }
  
  return(list(vectorized_Sigmas = vectorized_Sigmas, group_labels = labels))
}

#' Load data from Johnson et al. (2019)
#' 
#' @param visualize optional parameter specifying the quick visualization of the data to use;
#' "ribbon" fits and visualizes the fido model, "bars" plots a stacked bar chart
#' @param GP if TRUE, fit a Gaussian process model as opposed to a dynamic linear model
#' @return correlation matrix (hosts x taxon-taxon associations) or NULL
#' @export
load_Johnson_2D_summary <- function(visualize = NULL, GP = FALSE) {
  data_file <- file.path("input", "fit_johnson.rds")
  if(!file.exists(data_file) | !is.null(visualize)) {
    # Parse read count table
    counts <- read.table(file.path("input", "johnson2019", "taxonomy_counts_s.txt"),
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # Parse read sample ID to subject ID mapping
    mapping <- read.table(file.path("input", "johnson2019", "SampleID_map.txt"),
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Chop out taxonomy
    counts <- counts[,2:ncol(counts)]
    
    # Filter to some minimum relative average abundance
    retain_idx <- filter_taxa_2D_summary(counts)
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
    
    # # Testing
    # for(i in 1:length(host_columns)) {
    #   limit <- min(5, length(host_columns[[i]]))
    #   host_columns[[i]] <- host_columns[[i]][1:limit]
    #   host_dates[[i]] <- host_dates[[i]][1:limit]
    # }
    
    vectorized_Sigmas <- fit_or_visualize_sub(counts, host_columns, host_dates, visualize, data_file, GP)
  } else {
    vectorized_Sigmas <- readRDS(data_file)
  }
  return(vectorized_Sigmas)
}

#' Load data from Dethlefsen & Relman (2011)
#' 
#' @param visualize optional parameter specifying the quick visualization of the data to use;
#' "ribbon" fits and visualizes the fido model, "bars" plots a stacked bar chart
#' @param GP if TRUE, fit a Gaussian process model as opposed to a dynamic linear model
#' @return correlation matrix (hosts x taxon-taxon associations) or NULL
#' @import lubridate
#' @export
load_DethlefsenRelman_2D_summary <- function(visualize = NULL, GP = FALSE) {
  data_file <- file.path("input", "fit_dethlefsen.rds")
  if(!file.exists(data_file) | !is.null(visualize)) {
    # Read count table
    counts <- read.table(file.path("input", "dethlefsen_relman2011", "sd01.txt"),
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    counts <- counts[,2:163]
    
    # Filter to some minimum relative average abundance
    retain_idx <- filter_taxa_2D_summary(counts)
    counts <- counts[retain_idx,]
    
    # Read the sampling schedule metadata
    metadata <- read.table(file.path("input", "dethlefsen_relman2011", "sampling_schedule.txt"),
                           header = FALSE, stringsAsFactors = FALSE)
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
    
    # # Testing
    # for(i in 1:length(host_columns)) {
    #   limit <- min(5, length(host_columns[[i]]))
    #   host_columns[[i]] <- host_columns[[i]][1:limit]
    #   host_dates[[i]] <- host_dates[[i]][1:limit]
    # }
    
    vectorized_Sigmas <- fit_or_visualize_sub(counts, host_columns, host_dates, visualize, data_file, GP)
  } else {
    vectorized_Sigmas <- readRDS(data_file)
  }
  return(vectorized_Sigmas)
}

#' Load data from Caporaso et al. (2011)
#' 
#' @param visualize optional parameter specifying the quick visualization of the data to use;
#' "ribbon" fits and visualizes the fido model, "bars" plots a stacked bar chart
#' @param GP if TRUE, fit a Gaussian process model as opposed to a dynamic linear model
#' @return correlation matrix (hosts x taxon-taxon associations) or NULL
#' @export
load_Caporaso_2D_summary <- function(visualize = NULL, GP = FALSE) {
  data_file <- file.path("input", "fit_caporaso.rds")
  if(!file.exists(data_file) | !is.null(visualize)) {
    # Read count tables; OTUs all agree -- already checked
    counts_1 <- read.table(file.path("input", "caporaso2011", "F4_feces_L6.txt"),
                           header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    counts_1 <- counts_1[,2:ncol(counts_1)]
    
    counts_2 <- read.table(file.path("input", "caporaso2011", "M3_feces_L6.txt"),
                           header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    counts_2 <- counts_2[,2:ncol(counts_2)]
    
    # First row are the sample day identifiers
    counts <- cbind(counts_1, counts_2)
    
    # Filter to some minimum relative average abundance
    retain_idx <- filter_taxa_2D_summary(counts[2:nrow(counts),])
    counts <- counts[c(TRUE, retain_idx),]
    
    host_columns <- list(1:ncol(counts_1), (ncol(counts_1)+1):ncol(counts)) # individuals 1, 2
    host_dates <- list()
    for(i in 1:2) {
      host_dates[[i]] <- unname(unlist(counts[1,host_columns[[i]]]))
      host_dates[[i]] <- host_dates[[i]] - min(host_dates[[i]]) # preserve 0 start
    }
    
    # # Testing
    # for(i in 1:length(host_columns)) {
    #   limit <- min(5, length(host_columns[[i]]))
    #   host_columns[[i]] <- host_columns[[i]][1:limit]
    #   host_dates[[i]] <- host_dates[[i]][1:limit]
    # }
    
    vectorized_Sigmas <- fit_or_visualize_sub(counts, host_columns, host_dates, visualize, data_file, GP)
  } else {
    vectorized_Sigmas <- readRDS(data_file)
  }
  return(vectorized_Sigmas)
}

#' Load data from David et al. (2014)
#' 
#' @param visualize optional parameter specifying the quick visualization of the data to use;
#' "ribbon" fits and visualizes the fido model, "bars" plots a stacked bar chart
#' @param GP if TRUE, fit a Gaussian process model as opposed to a dynamic linear model
#' @return correlation matrix (hosts x taxon-taxon associations) or NULL
#' @import stringr
#' @export
load_David_2D_summary <- function(visualize = NULL, GP = FALSE) {
  data_file <- file.path("input", "fit_caporaso.rds")
  if(!file.exists(data_file) | !is.null(visualize)) {
    counts <- read.table(file.path("input", "david2014", "otu.table.ggref"),
                         header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    metadata <- read.table(file.path("input", "david2014", "13059_2013_3286_MOESM18_ESM.csv"),
                           sep = ",", header = TRUE)
    
    sample_labels <- metadata$X
    subj_labels <- metadata$AGE # subject A is 26, subject B is 36
    day_labels <- metadata$COLLECTION_DAY
    subj_labels[subj_labels == 26] <- "A"
    subj_labels[subj_labels == 36] <- "B"
    # Strip character from the count column names
    sample_labels <- sapply(sample_labels, function(x) {
      str_replace(x, "\\.\\d+", "")
    })
    sample_labels <- unname(sample_labels)
    
    # Remove saliva samples
    keep_idx <- which(!sapply(sample_labels, function(x) {
      str_detect(x, "Saliva")
    }))
    subj_labels <- subj_labels[keep_idx]
    sample_labels <- sample_labels[keep_idx]
    day_labels <- day_labels[keep_idx]
    
    # Include columns in counts that are in sample_labels
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
    # Filter to some minimum relative average abundance
    retain_idx <- filter_taxa_2D_summary(counts)
    counts <- counts[retain_idx,]
    counts <- as.matrix(counts)
    
    host_columns <- list(1:length(subj_A_days), (length(subj_A_days) + 1):ncol(counts)) # individuals 1, 2
    host_dates <- list(subj_A_days, subj_B_days)
    
    # # Testing
    # for(i in 1:length(host_columns)) {
    #   limit <- min(5, length(host_columns[[i]]))
    #   host_columns[[i]] <- host_columns[[i]][1:limit]
    #   host_dates[[i]] <- host_dates[[i]][1:limit]
    # }
    
    vectorized_Sigmas <- fit_or_visualize_sub(counts, host_columns, host_dates, visualize, data_file, GP)
  } else {
    vectorized_Sigmas <- readRDS(data_file)
  }
  return(vectorized_Sigmas)
}

#' Load data from DIABIMMUNE study
#' 
#' @param visualize optional parameter specifying the quick visualization of the data to use;
#' "ribbon" fits and visualizes the fido model, "bars" plots a stacked bar chart
#' @param GP if TRUE, fit a Gaussian process model as opposed to a dynamic linear model
#' @return correlation matrix (hosts x taxon-taxon associations) or NULL
#' @export
load_DIABIMMUNE_2D_summary <- function(visualize = NULL, GP = FALSE) {
  data_file <- file.path("input", "fit_DIABIMMUNE.rds")
  if(!file.exists(data_file) | !is.null(visualize)) {
    load(file.path("input", "DIABIMMUNE", "diabimmune_karelia_16s_data.rdata"))
    load(file.path("input", "DIABIMMUNE", "DIABIMMUNE_Karelia_metadata.RData"))
    counts <- t(data_16s)
    rm(data_16s)
    
    # these aren't counts but relative abundances
    # let's dream and scale them into counts
    for(i in 1:ncol(counts)) {
      counts[,i] <- rmultinom(1, size = rpois(1, 10000), prob = counts[,i])
    }
    retain_idx <- filter_taxa_2D_summary(counts)
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
      sample_days <- sapply(host_dates[[i]], function(x) x - baseline_day)
      host_dates[[i]] <- sample_days
    }
    
    # (Laboriously) fix the date ordering
    for(host in names(host_columns)) {
      reorder <- order(host_dates[[host]])
      host_dates[[host]] <- host_dates[[host]][reorder]
      host_columns[[host]] <- host_columns[[host]][reorder]
    }
    
    # # Testing
    # for(i in 1:length(host_columns)) {
    #   limit <- min(5, length(host_columns[[i]]))
    #   host_columns[[i]] <- host_columns[[i]][1:limit]
    #   host_dates[[i]] <- host_dates[[i]][1:limit]
    # }
  
    vectorized_Sigmas <- fit_or_visualize_sub(counts, host_columns, host_dates, visualize, data_file, GP)
  } else {
    vectorized_Sigmas <- readRDS(data_file)
  }
  return(vectorized_Sigmas)
}

#' Load data from Grossart German lakes study
#' 
#' @param visualize optional parameter specifying the quick visualization of the data to use;
#' "ribbon" fits and visualizes the fido model, "bars" plots a stacked bar chart
#' @param GP if TRUE, fit a Gaussian process model as opposed to a dynamic linear model
#' @return correlation matrix (hosts x taxon-taxon associations) or NULL
#' @import phyloseq
#' @export
load_Grossart_2D_summary <- function(visualize = NULL, GP = FALSE) {
  data_file <- file.path("input", "fit_grossart.rds")
  if(!file.exists(data_file) | !is.null(visualize)) {
    counts <- readRDS(file.path("input", "grossart_lakes", "data.rds"))
    counts <- otu_table(counts)@.Data
    retain_idx <- filter_taxa_2D_summary(counts)
    counts <- counts[retain_idx,]
    
    metadata <- read.table("input/grossart_lakes/945_20190102-074005.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    # Check out the sample numbers per unique combo of depth and filter poresize for each lake
    # Below, we'll hand-pick similar conditions with large sample numbers!
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
      days_vector <- unname(sapply(timestamps, function(x) difftime(as.Date(x), as.Date(baseline_date), units = "days")))
      host_columns[[i]] <- which(colnames(counts) %in% sample_names)
      host_dates[[i]] <- days_vector
    }
    
    # # Testing
    # for(i in 1:length(host_columns)) {
    #   limit <- min(5, length(host_columns[[i]]))
    #   host_columns[[i]] <- host_columns[[i]][1:limit]
    #   host_dates[[i]] <- host_dates[[i]][1:limit]
    # }

    vectorized_Sigmas <- fit_or_visualize_sub(counts, host_columns, host_dates, visualize, data_file, GP)
  } else {
    vectorized_Sigmas <- readRDS(data_file)
  }
  return(vectorized_Sigmas)
}

#' Load epilimnion (shallow water) data from McMahon Wisconsin, US lakes study
#' 
#' @param visualize optional parameter specifying the quick visualization of the data to use;
#' "ribbon" fits and visualizes the fido model, "bars" plots a stacked bar chart
#' @param GP if TRUE, fit a Gaussian process model as opposed to a dynamic linear model
#' @return correlation matrix (hosts x taxon-taxon associations) or NULL
#' @import phyloseq
#' @import stringr
#' @export
load_McMahon_shallow_2D_summary <- function(visualize = NULL, GP = FALSE) {
  data_file <- file.path("input", "fit_mcmahon_shallow.rds")
  if(!file.exists(data_file) | !is.null(visualize)) {
    counts <- readRDS(file.path("input", "mcmahon_lakes", "data.rds"))
    counts <- otu_table(counts)@.Data
    retain_idx <- filter_taxa_2D_summary(counts)
    counts <- counts[retain_idx,]
    
    metadata <- read.table(file.path("input", "mcmahon_lakes", "1288_20180418-110149.txt"),
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    metadata <- parse_McMahon_metadata(metadata)
    
    site_obj <- parse_McMahon_sites(counts, metadata, layer = "E")
    host_columns <- site_obj$host_columns
    host_dates <- site_obj$host_dates
    # There are duplicate samples on some days but these look legit (?)
    
    # # Testing
    # for(i in 1:length(host_columns)) {
    #   limit <- min(5, length(host_columns[[i]]))
    #   host_columns[[i]] <- host_columns[[i]][1:limit]
    #   host_dates[[i]] <- host_dates[[i]][1:limit]
    # }
    
    vectorized_Sigmas <- fit_or_visualize_sub(counts, host_columns, host_dates, visualize, data_file, GP)
  } else {
    vectorized_Sigmas <- readRDS(data_file)
  }
  return(vectorized_Sigmas)
}


#' Load hypolimnion (deep water) data from McMahon Wisconsin, US lakes study
#' 
#' @param visualize optional parameter specifying the quick visualization of the data to use;
#' "ribbon" fits and visualizes the fido model, "bars" plots a stacked bar chart
#' @param GP if TRUE, fit a Gaussian process model as opposed to a dynamic linear model
#' @return correlation matrix (hosts x taxon-taxon associations) or NULL
#' @import phyloseq
#' @import stringr
#' @export
load_McMahon_deep_2D_summary <- function(visualize = NULL, GP = FALSE) {
  data_file <- file.path("input", "fit_mcmahon_deep.rds")
  if(!file.exists(data_file) | !is.null(visualize)) {
    counts <- readRDS(file.path("input", "mcmahon_lakes", "data.rds"))
    counts <- otu_table(counts)@.Data
    retain_idx <- filter_taxa_2D_summary(counts)
    counts <- counts[retain_idx,]
    
    metadata <- read.table(file.path("input", "mcmahon_lakes", "1288_20180418-110149.txt"),
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    
    metadata <- parse_McMahon_metadata(metadata)
    
    site_obj <- parse_McMahon_sites(counts, metadata, layer = "H")
    host_columns <- site_obj$host_columns
    host_dates <- site_obj$host_dates
    
    # # Testing
    # for(i in 1:length(host_columns)) {
    #   limit <- min(5, length(host_columns[[i]]))
    #   host_columns[[i]] <- host_columns[[i]][1:limit]
    #   host_dates[[i]] <- host_dates[[i]][1:limit]
    # }
    
    vectorized_Sigmas <- fit_or_visualize_sub(counts, host_columns, host_dates, visualize, data_file, GP)
  } else {
    vectorized_Sigmas <- readRDS(data_file)
  }
  return(vectorized_Sigmas)
}

parse_McMahon_metadata <- function(metadata) {
  metadata$lake <- NA
  metadata$layer <- NA
  # Append more useful columns to the metadata
  for(i in 1:nrow(metadata)) {
    matches <- str_match_all(metadata$description[i], regex("^freshwater metagenome (\\D+)(\\d+)(\\D+)(\\d+)"))
    lake_name <- matches[[1]][1,2]
    if(lake_name != "NSb") {
      layer <- str_match_all(lake_name, regex("[H|E]$"))
      if(nrow(layer[[1]]) > 0) {
        if(layer[[1]][1,1] == "H") {
          # hypolimnion; deep water
          metadata$lake[i] <- substr(lake_name, 1, str_length(lake_name)-1)
          metadata$layer[i] <- "H"
        } else {
          # epilimnion; shallow water
          metadata$lake[i] <- substr(lake_name, 1, str_length(lake_name)-1)
          metadata$layer[i] <- "E"
        }
      } else {
        metadata$lake[i] <- lake_name
      }
    }
  }
  return(metadata)
}

parse_McMahon_sites <- function(counts, metadata, layer) {
  host_columns <- list()
  host_dates <- list()
  for(i in 1:ncol(counts)) {
    sample_id <- colnames(counts)[i]
    metadata_idx <- which(metadata$sample_name == sample_id)
    lake_name <- metadata$lake[metadata_idx]
    layer_label <- metadata$layer[metadata_idx]
    if(!is.na(layer_label) & layer_label == layer & !is.na(lake_name)) { # "NSb", the singleton
      sample_date <- as.Date(metadata$collection_timestamp[metadata_idx])
      if(lake_name %in% names(host_columns)) {
        host_columns[[lake_name]] <- c(host_columns[[lake_name]], i)
        host_dates[[lake_name]] <- c(host_dates[[lake_name]], sample_date)
      } else {
        host_columns[[lake_name]] <- c(i)
        host_dates[[lake_name]] <- sample_date
      }
    }
  }
  # Convert dates to days
  for(i in 1:length(host_dates)) {
    baseline_date <- min(host_dates[[i]])
    sample_days <- sapply(host_dates[[i]], function(x) difftime(x, baseline_date, units = "days"))
    host_dates[[i]] <- sample_days
  }
  # (Laboriously) fix the date ordering
  for(i in 1:length(host_columns)) {
    reorder <- order(host_dates[[i]])
    host_dates[[i]] <- host_dates[[i]][reorder]
    host_columns[[i]] <- host_columns[[i]][reorder]
  }
  return(list(host_columns = host_columns, host_dates = host_dates))
}

fit_or_visualize_sub <- function(counts, host_columns, host_dates, visualize, data_file, GP = FALSE) {
  if(!is.null(visualize)) {
    subject_idx <- 1
    # counts <- counts[,c(1, host_columns[[subject_idx]])]
    host_columns <- list(host_columns[[subject_idx]])
    host_dates <- list(host_dates[[subject_idx]])
    if(visualize == "ribbon") {
      fit <- fit_model(counts, host_columns, host_dates, depth = 100, GP = GP)
      visualize_fit(fit, counts[,host_columns[[1]]], host_dates[[1]], GP = GP)
    } else {
      plot_bars(counts[,host_columns[[1]]])
    }
    return(NULL)
  } else {
    vectorized_Sigmas <- fit_model(counts, host_columns, host_dates, GP = GP)
  }
  vectorized_Sigmas <- vectorized_Sigmas[,,1]
  saveRDS(vectorized_Sigmas, file = data_file)
  return(vectorized_Sigmas)
}

#' @import fido
fit_model <- function(counts, host_columns, host_dates, depth = 1, GP = FALSE) {
  if(depth < 1) {
    depth <- 1
  }
  D <- nrow(counts)
  n_interactions <- n_pairwise_combos(D)
  vectorized_Sigmas <- array(NA, dim = c(length(host_columns), n_interactions, depth))
  for(subject in 1:length(host_columns)) {
    subject_samples <- host_columns[[subject]]
    subject_dates <- host_dates[[subject]] + 1
    if(length(subject_samples) > 0) {
      subject_counts <- counts[,subject_samples]
      Y <- as.matrix(subject_counts)
      
      # Check that the dates are in chronological order
      # (The DLM implementation assumes this)
      # reorder <- order(subject_dates)
      # Y <- Y[,reorder]
      # subject_dates <- subject_dates[reorder]
      
      alr_ys <- driver::alr((t(Y) + 0.5))
      alr_means <- colMeans(alr_ys)
      taxa_covariance <- get_Xi(D, total_variance = 1)
      X <- matrix(subject_dates, 1, length(subject_dates))
      
      if(depth == 1) {
        n_samples <- 0
        ret_mean <- TRUE
      } else {
        n_samples <- depth
        ret_mean <- FALSE
      }
      
      if(GP) {
        Theta <- function(X) matrix(alr_means, D-1, ncol(X))
        rho <- calc_se_decay(min_correlation = 0.1, days_to_baseline = 7)
        Gamma <- function(X) {
          SE(X, sigma = 1, rho = rho, jitter = 1e-08)
        }
        fit <- basset(Y, X, taxa_covariance$upsilon, Theta, Gamma, taxa_covariance$Xi,
                            n_samples = n_samples, ret_mean = ret_mean)
      } else {
        # DLM; set params as in fit_DLM()
        T <- max(subject_dates)
        F <- matrix(1, 1, T)
        Q <- nrow(F)
        W <- diag(Q)
        W <- W/nrow(W)
        G <- diag(Q)
        C0 <- W
        M0 <- matrix(0, Q, D-1)
        M0[1,] <- alr_means
        fit <- labraduck(Y = Y, upsilon = taxa_covariance$upsilon, Xi = taxa_covariance$Xi,
                         F = F, G = G, W = W, M0 = M0, C0 = C0,
                         observations = X, gamma_scale = 1,
                         W_scale = 1, apply_smoother = (depth > 1),
                         n_samples = n_samples, ret_mean = ret_mean,
                         pars = c("Eta", "Sigma", "Thetas_filtered", "Thetas_smoothed", "Eta_DLM"))
      }
      
      if(depth > 1) {
        return(fit)
      }
      fit.clr <- to_clr(fit)
      Sigmas <- fit.clr$Sigma
      for(k in 1:depth) {
        temp <- cov2cor(Sigmas[,,k])
        temp <- c(temp[upper.tri(temp, diag = FALSE)])
        vectorized_Sigmas[subject,,k] <- temp
      }
    }
  }
  return(vectorized_Sigmas)
}

#' @import driver
#' @import dplyr
#' @import ggplot2
visualize_fit <- function(fit, subject_counts, subject_dates, GP = FALSE) {
  if(fit$iter > 1) {
    subject_dates <- subject_dates + 1
    D <- nrow(subject_counts)
    coord <- sample(1:D)[1] # Randomly pick a logratio taxon to plot
    if(GP) {
      # GP -- renders CLR trajectory (kind of arbitrarily); should probably change reference points to
      # inferred Etas, not directly transformed observations (TBD)
      lr_ys <- clr(t(fit$Y) + 0.5)
      lr_coord <- data.frame(x = subject_dates, y = lr_ys[,coord])
      Eta <- predict(fit, min(subject_dates):max(subject_dates), response = "Eta", iter = fit$iter)
      Eta <- alrInv_array(Eta, fit$D, 1)
      Eta <- clr_array(Eta, 1)
      posterior_samples <- gather_array(Eta[coord,,], "logratio_value", "observation", "sample_number")
    } else {
      # DLM -- renders ALR trajectory
      lr_ys <- t(apply(fit$Eta, c(1,2), mean))
      lr_coord <- data.frame(x = subject_dates, y = lr_ys[,coord])
      F <- fit$F
      Y <- fit$Y
      T <- max(subject_dates)
      N <- length(subject_dates)
      n_samples <- dim(fit$Eta)[3]
      Ft <- t(F)
      Q <- nrow(F)
      Theta_samples <- matrix(NA, T, n_samples)
      Eta_samples <- matrix(NA, T, n_samples)
      # Iterate through the posterior samples
      # The fitted model contains the interpolated time points but these are 3D
      # passed back as 2D matrices and need to be reshaped
      for(k in 1:n_samples) {
        ThetasS_1T <- fit$Thetas_smoothed[,k]
        dim(ThetasS_1T) <- c(Q, D-1, T)
        Theta_samples[,k] <- ThetasS_1T[,coord,]
        EtasS_1T <- fit$Eta_DLM[,k]
        dim(EtasS_1T) <- c(Q, D-1, T)
        Eta_samples[,k] <- EtasS_1T[,coord,]
      }
      posterior_samples <- gather_array(Eta_samples, "logratio_value", "observation", "sample_number")
    }
    # Get quantiles
    post_quantiles <- posterior_samples %>%
      group_by(observation) %>%
      summarise(p2.5 = quantile(logratio_value, prob=0.025),
                p25 = quantile(logratio_value, prob=0.25),
                mean = mean(logratio_value),
                p75 = quantile(logratio_value, prob=0.75),
                p97.5 = quantile(logratio_value, prob=0.975)) %>%
      ungroup()
    p <- ggplot(post_quantiles, aes(x=observation, y=mean)) +
      geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
      geom_ribbon(aes(ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
      geom_line(color="blue") +
      geom_point(data = lr_coord, aes(x = x, y = y), alpha=0.5) +
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle=45)) +
      ylab("LR coord")
    show(p)
  }
}

#' Plot relative abundances as stacked bars
#' 
#' @param counts count table (taxa x samples)
#' @param palette optional color palette (as a vector of hex color strings)
#' @return NULL
#' @import driver
#' @import ggplot2
#' @export
plot_bars <- function(counts, palette = NULL) {
  df <- gather_array(counts, "value", "taxon", "time")
  df$taxon <- factor(df$taxon)
  if(is.null(palette)) {
    palette <- generate_highcontrast_palette(length(levels(df$taxon)))
  }
  p <- ggplot(df) + 
    geom_bar(aes(x = time, y = value, fill = taxon), position = "fill", stat = "identity") +
    scale_fill_manual(values = palette) +
    theme(legend.position="none")
  show(p)
}





