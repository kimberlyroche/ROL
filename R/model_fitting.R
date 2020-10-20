#' Periodic kernel
#' 
#' @param X covariate (dimension Q x N; i.e., covariates x samples) 
#' @param sigma scalar parameter 
#' @param rho scalar bandwidth parameter
#' @param period period length (in days)
#' @param jitter small scalar to add to off-diagonal of gram matrix 
#'   (for numerical underflow issues)
#' @return Gram Matrix (N x N) (e.g., the Kernel evaluated at 
#' each pair of points)
#' @export
PER <- function(X, sigma=1, rho=1, period=24, jitter=0){
  dist <- as.matrix(dist(t(X)))
  G <- sigma^2 * exp(-2*(sin(pi*dist/period)^2)/(rho^2)) + jitter*diag(ncol(dist))
}

#' Trivially expands the dimensions of the posterior estimates objects in a MAP estimate basset result
#' 
#' @param fit a bassetfit object
#' @return NULL
#' @export
#' @examples
#' fit <- fix_MAP_dims(fit)
fix_MAP_dims <- function(fit) {
  if(length(dim(fit$fit$Eta)) == 2) {
    dim(fit$fit$Eta) <- c(dim(fit$fit$Eta),1)
    dim(fit$fit$Lambda) <- c(dim(fit$fit$Lambda),1)
    dim(fit$fit$Sigma) <- c(dim(fit$fit$Sigma),1)
  }
  return(fit)
}

#' Pulls a list of fitted basset models from the designated output directory
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param DLM if TRUE, looks for DLM model fits instead of GP model fits
#' @param MAP use MAP estimate model output instead of full posterior output
#' @return NULL
#' @export
#' @examples
#' model_list <- get_fitted_model_list(tax_level = "ASV", DLM = TRUE, MAP = FALSE)
get_fitted_model_list <- function(tax_level = "ASV", DLM = FALSE, MAP = FALSE) {
  if(DLM) {
    pattern_str <- "*_labraduckfit.rds"
    regexpr_str <- "_labraduckfit.rds"
  } else {
    pattern_str <- "*_bassetfit.rds"
    regexpr_str <- "_bassetfit.rds"
  }
  if(MAP) {
    level_dir <- file.path("output","model_fits",paste0(tax_level,"_MAP"))
  } else {
    level_dir <- file.path("output","model_fits",tax_level)
  }
  fitted_models <- list.files(path=level_dir, pattern=pattern_str, full.names=TRUE, recursive=FALSE)
  # pull out some useful summary information: a list of fitted hosts, models, and dimensions  
  hosts <- as.vector(sapply(fitted_models, function(x) { idx <- regexpr(regexpr_str, x); return(substr(x, idx-3, idx-1)) } ))
  if(MAP) {
    path <- file.path("output","model_fits",paste0(tax_level,"_MAP/",hosts[1],regexpr_str))
  } else {
    path <- file.path("output","model_fits",paste0(tax_level,"/",hosts[1],regexpr_str))
  }
  fit <- read_file(path)
  return(list(hosts=hosts,
              pattern_str=pattern_str,
              regexpr_str=regexpr_str,
              model_list=fitted_models,
              D=fit$fit$D,
              n_samples=fit$fit$iter))
}

#' Set up a basic ALR prior
#' 
#' @param D number of features including reference (where the ALR will represent D-1)
#' @param total_variance scale of the log variance
#' @return list containing inverse Wishart parameters degrees of freedom and scale matrix
#' @export
#' @examples
#' params <- get_Xi(D=100, log_var_scale=1)
#' Sigma <- matrixsampling::rinvwishart(1, params$upsilon, params$Xi)[,,1]
get_Xi <- function(D, total_variance=1) {
  upsilon <- D-1+10 # specify low certainty/concentration
  GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference
  Xi <- GG%*%(diag(D)*total_variance)%*%t(GG) # take diag as covariance over log abundances
  Xi <- Xi*(upsilon-D-1)
  return(list(upsilon=upsilon, Xi=Xi))
}

#' Generate empirical estimates for variance components to use in the Gaussian process
#' 
#' @param data a phyloseq object
#' @return list containing median log abundance taxon and mean empirical total variance
#' across hosts using this ALR representation
#' @import phyloseq
#' @import driver
#' @export
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

#' Build the design matrix
#' 
#' @param metadata metadata data.frame associated with the data phyloseq object
#' @param metadata.diet metadata data.frame parsed from separate diet and climate data
#' @return matrix with observation time (days), diet PCs, and climate (rainfall and max temp) covariates
#' @export
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

#' Build the design matrix for prediction by interpolating within the 
#' 
#' @param X existing design matrix
#' @return matrix with observation time (days), diet PCs, and climate (rainfall and max temp) covariates interpolated from
#' earliest to latest observation date
#' @export
build_design_matrix_predict <- function(X) {
  last_time_point <- max(X[1,])
  X_predict <- t(1:last_time_point)
  n_features <- nrow(X)-1 # these are the features to (linearly) interpolate
  if(n_features > 0) {
    for(f in 1:n_features) {
      X_predict <- rbind(X_predict, t(approx(x = X[1,], y = X[f+1,], xout = 1:last_time_point, method = "linear")$y))
    }
  }
  return(X_predict)
}

#' Define bandwidth of squared exponential kernel
#' 
#' @param min_correlation minimum correlation to assume between (within-host) samples
#' @param days_to_baseline days at which squared exponential kernel decays to baseline correlation of ~0.1
#' @return bandwidth parameter for squared exponential kernel
#' @export
calc_se_decay <- function(min_correlation = 0.1, days_to_baseline = 90) {
  # back-calculate the squared exponential bandwidth parameter by finding a bandwidth that gives
  # a desired minimum correlation at the number of days specified by days_to_baseline  
  rho <- sqrt(-days_to_baseline^2/(2*log(min_correlation))) # back calculate the decay
  return(rho)
}

#' Define a kernel (function) over samples
#' 
#' @param kernel_scale total variance for the composite kernel
#' @param proportions proportion variance to attribute to each of 3 kernels (see details)
#' @param rho bandwidth for SE kernel
#' @details Composite kernel is built from (1) squared exponential kernel (base autocorrelation component),
#' (2) seasonal kernel (periodic), and (3) diet and climate component (another squared exponential)
#' @return list containing kernel function and bandwidth parameter
#' @import fido
#' @export
get_Gamma <- function(kernel_scale, proportions, min_correlation = 0.1, days_to_baseline = 90) {
  rho <- calc_se_decay(min_correlation = min_correlation, days_to_baseline = days_to_baseline)
  # back-calculate the squared exponential bandwidth parameter by finding a bandwidth that gives
  # a desired minimum correlation at the number of days specified by SE_days_to_baseline
  Gamma <- function(X) {
    jitter <- 1e-08
    proportions <- abs(proportions)/sum(abs(proportions)) # just in case we pass something that isn't a composition
    part.1 <- kernel_scale * proportions[1]
    part.2 <- kernel_scale * proportions[2]
    part.3 <- kernel_scale * proportions[3]
    SE(X[1,,drop=F], sigma = sqrt(part.1), rho = rho, jitter = jitter) +
      PER(X[1,,drop=F], sigma = sqrt(part.2), rho = 1, period = 365, jitter = jitter) +
      SE(X[2:6,,drop=F], sigma = sqrt(part.3), rho = 1, jitter = jitter) # a few diet PCs
  }
  return(list(rho = rho, Gamma = Gamma))
}

#' Utility function to calculate RMSE
#'
#' @param y1 count vector
#' @param y2 count vector
#' @param pc optional pseudocount parameter
#' @return root mean squared error of log counts
#' @export
log_rmse <- function(y1, y2, pc = 0.5) {
  sqrt(sum(sapply(1:length(y1), function(yy) {
    (log(y1[yy] + pc) - log(y2[yy] + pc))^2
  }))/length(y1))
}

#' Utility function to calculate fold change (as error)
#'
#' @param y1 count vector
#' @param y2 count vector
#' @param pc optional pseudocount parameter
#' @return average fold difference between vectors
#' @export
fold_error <- function(y1, y2, pc = 0.5) {
  mean(sapply(1:length(y1), function(yy) {
    a <- max(y1[yy] + pc, y2[yy] + pc)
    b <- min(y1[yy] + pc, y2[yy] + pc)
    a/b
  })) - 1
}

#' Fit a Gaussian process to a single host series using basset
#' 
#' @param data a phyloseq object
#' @param host host short name (e.g. ACA)
#' @param taxa_covariance list of prior covariance parameters over taxa
#' @param sample_covariance kernel function
#' @param rho bandwidth for SE kernel
#' @param tax_level taxonomic level at which to agglomerate data
#' @param alr_ref index of reference ALR coordinate
#' @param n_samples number of posterior samples to draw
#' @param MAP compute MAP estimate only (as single posterior sample)
#' @param holdout_proportion if non-zero, proportion of host's sample to use as a test set
#' @param return_model if TRUE, returns the fitted model instead of saving it
#' @param scramble if TRUE, host samples are scrambled in time; this serve the diagnostic purpose of 
#'        identifying whether some of the correlation we see between individual hosts' dynamics
#'        are the result of systematic differences in abundances (etc.)
#' @details Fitted model and metadata saved to designated model output directory.
#' @return NULL
#' @import phyloseq
#' @import fido
#' @export
#' @examples
#' tax_level <- "ASV"
#' data <- load_data(tax_level = tax_level)
#' params <- formalize_parameters(data)
#' sample_covariance <- get_Gamma(kernel_scale = 2, proportions = c(1, 0, 0), min_correlation = 0.1, days_to_baseline = 90)
#' taxa_covariance <- get_Xi(phyloseq::ntaxa(data), total_variance = 1)
#' fit_GP(data, host = "GAB", taxa_covariance = taxa_covariance, sample_covariance = sample_covariance, tax_level = tax_level, alr_ref = params$alr_ref, MAP = TRUE)
fit_GP <- function(data, host, taxa_covariance, sample_covariance, tax_level = "ASV", alr_ref = NULL, n_samples = 100, MAP = FALSE, holdout_proportion = 0, return_model = FALSE, scramble = FALSE) {
  if(MAP) {
    cat(paste0("Fitting fido::basset model (MAP) to host ",host,"\n"))
    if(holdout_proportion > 0) {
      stop("Leave-one-out predictions only work for full posterior estimates (currently).\n")
    }
  } else {
    cat(paste0("Fitting fido::basset model to host ",host,"\n"))
  }
  
  # global assign is a hack seemingly necessary for this phyloseq::subset_samples function call
  host <<- host
  host_data <- subset_samples(data, sname == host)
  
  # encode observations as differences from baseline in units of days
  host_metadata <- sample_data(host_data)

  # read diet and climate covariate data
  data.diet <- readRDS("input/ps_w_covs.RDS")
  data.name_mapping <- read.csv("input/host_subject_id_to_sname_key.csv")
  data.name_mapping <- unique(data.name_mapping[,c("sname","host_subject_id2")])
  host.num <<- as.character(data.name_mapping[data.name_mapping$sname == host,]$host_subject_id2)
  data.diet <- subset_samples(data.diet, host == host.num)
  metadata.diet <- sample_data(data.diet)

  # pull out the count table
  Y <- otu_table(host_data)@.Data
  X <- build_design_matrix(host_metadata, metadata.diet) # covariates x samples
  Y <- t(Y)                                              # taxa x samples
  
  # get dimensions
  D <- nrow(Y)
  N <- ncol(Y)

  # scramble, if required
  if(scramble) {
    cat("Scrambling host observations...\n")
    for(i in 1:D) {
      per_taxon_col_assignments <- sample(1:N)
      Y[i,] <- Y[i,per_taxon_col_assignments]
    }
  }
  
  # fido uses the D^th element as the ALR reference by default
  # if we'd like to use a different reference, do some row shuffling in Y to put the reference at the end
  if(!is.null(alr_ref)) {
    Y <- Y[c(setdiff(1:D,alr_ref),alr_ref),]
  }

  # strip off sequence variant labels
  colnames(Y) <- NULL
  rownames(Y) <- NULL

  if(holdout_proportion > 0) {
    holdout_n <- round(N*holdout_proportion)
    holdout_idx <- sample(1:N)[1:holdout_n]
    X_test <- X[,holdout_idx]
    X <- X[,setdiff(1:N, holdout_idx)]
    Y_test <- Y[,holdout_idx]
    Y <- Y[,setdiff(1:N, holdout_idx)]
  }

  # define the prior over baselines
  alr_ys <- driver::alr((t(Y) + 0.5))
  alr_means <- colMeans(alr_ys)
  Theta <- function(X) matrix(alr_means, D-1, ncol(X))
  
  if(MAP) {
    n_samples <- 0
  }
  fit <- fido::basset(Y, X, taxa_covariance$upsilon, Theta, sample_covariance$Gamma, taxa_covariance$Xi,
                       n_samples = n_samples, ret_mean = MAP, 
                       b2 = 0.98, step_size = 0.004, eps_f = 1e-11, eps_g = 1e-05,
                       max_iter = 10000L, optim_method = "adam")

  if(MAP) {
    # fill out dimensions; some function (including fido predictive functions) expect arrays
    # not matrices here
    fit <- fix_MAP_dims(fit)
  }

  if(holdout_proportion > 0) {
    predicted <- predict(fit, X_test, response = "Y", iter = fit$iter) # taxa x holdout samples x posterior samples
    err1 <- mean(sapply(1:dim(predicted)[2], function(i) {
                   # for each held out sample
                   ref <- Y_test[,i]
                   sapply(1:dim(predicted)[3], function(j) {
                     # for each posterior sample
                     fold_error(ref, predicted[,i,j])
                   })
                 }))
    err2 <- mean(sapply(1:dim(predicted)[2], function(i) {
                   # for each held out sample
                   ref <- Y_test[,i]
                   sapply(1:dim(predicted)[3], function(j) {
                     # for each posterior sample
                     log_rmse(ref, predicted[,i,j])
                   })
                 }))
    return(list(fold_error = err1, log_rmse = err2))
  }
  
  # pack up results
  fit_obj <- list(Y = Y, alr_ys = alr_ys, alr_ref = alr_ref, X = X, fit = fit)
  
  # a hack: for later prediction, these will need to be available in the workspace for Gamma
  fit_obj$kernelparams$rho <- sample_covariance$rho

  if(return_model) {
    return(fit_obj)
  }

  # save results
  if(MAP) {
    save_dir <- check_output_dir(c("output","model_fits",paste0(tax_level,"_MAP")))
  } else {
    save_dir <- check_output_dir(c("output","model_fits",tax_level))
  }
  saveRDS(fit_obj, file.path(save_dir,paste0(host,"_bassetfit.rds")))
}

#' Fit a dynamic linear model to a single host series using labraduck
#' 
#' @param data a phyloseq object
#' @param host host short name (e.g. ACA)
#' @param taxa_covariance list of prior covariance parameters over taxa
#' @param var_scale combined scale of the total variance components Sigma and Gamma (default 1)
#' @param tax_level taxonomic level at which to agglomerate data
#' @param alr_ref index of reference ALR coordinate
#' @param n_samples number of posterior samples to draw
#' @param MAP compute MAP estimate only (as single posterior sample)
#' @details Fitted model and metadata saved to designated model output directory.
#' @return NULL
#' @import phyloseq
#' @import fido
#' @export
#' @examples
#' tax_level <- "ASV"
#' data <- load_data(tax_level = tax_level)
#' params <- formalize_parameters(data)
#' taxa_covariance <- get_Xi(phyloseq::ntaxa(data), total_variance = 1)
#' fit_DLM(data, host = "GAB", taxa_covariance = taxa_covariance, tax_level = tax_level, alr_ref = params$alr_ref, MAP = TRUE)
fit_DLM <- function(data, host, taxa_covariance, var_scale = 1, tax_level = "ASV", alr_ref = NULL, n_samples = 100, MAP = FALSE) {
  if(MAP) {
    cat(paste0("Fitting fido::labraduck model (MAP) to host ",host,"\n"))
  } else {
    cat(paste0("Fitting fido::labraduck model to host ",host,"\n"))
  }
  
  # global assign is a hack seemingly necessary for this phyloseq::subset_samples function call
  host <<- host
  host_data <- subset_samples(data, sname == host)
  
  # encode observations as differences from baseline in units of days
  host_metadata <- sample_data(host_data)

  # read diet and climate covariate data
  data.diet <- readRDS("input/ps_w_covs.RDS")
  data.name_mapping <- read.csv("input/host_subject_id_to_sname_key.csv")
  data.name_mapping <- unique(data.name_mapping[,c("sname","host_subject_id2")])
  host.num <<- as.character(data.name_mapping[data.name_mapping$sname == host,]$host_subject_id2)
  data.diet <- subset_samples(data.diet, host == host.num)
  metadata.diet <- sample_data(data.diet)

  days <- host_metadata$collection_date
  day0 <- min(days)
  days <- round(unname(sapply(days, function(x) difftime(x, day0, units = "days")))) + 1
  
  T <- max(days)
  # Build the pseudo-covariate matrix
  F <- matrix(0, 3, T)
  F[1,] <- 1
  for(i in 1:length(days)) {
    # F[2,data$days[i]] <- data$season[i] # season actually worsens the fit here
    F[2,days[i]] <- as.vector(scale(metadata.diet$PC1))[i]
    F[3,days[i]] <- as.vector(scale(metadata.diet$rain_monthly))[i]
  }

  # pull out the count table
  Y <- t(otu_table(host_data)@.Data) # taxa x samples
  # strip off sequence variant labels
  colnames(Y) <- NULL
  rownames(Y) <- NULL
  
  Q <- nrow(F)
  D <- nrow(Y)
  W <- diag(Q)
  # scale covariate-inclusive and covariate-exclusive models to have similar total variance
  W <- W/nrow(W)
  G <- diag(Q)

  # fido uses the D^th element as the ALR reference by default
  # if we'd like to use a different reference, do some row shuffling in Y to put the reference at the end
  if(!is.null(alr_ref)) {
    Y <- Y[c(setdiff(1:D,alr_ref),alr_ref),]
  }

  # define the prior over baselines
  C0 <- W
  alr_ys <- driver::alr((t(Y) + 0.5))
  alr_means <- colMeans(alr_ys)
  M0 <- matrix(0, Q, D-1)
  M0[1,] <- alr_means
  
  if(MAP) {
    n_samples <- 0
  }
  
  # I'm giving W about 1/2 the scale of gamma
  # My thinking here is that in the gLV models we've been simulating from, we've seen that most of
  #   the dynamism in the model results from /environmental interactions/, not innate volatility.
  #   In my mind, giving W a smaller scale than gamma is consistent with relatively more volatility
  #   presumed to result from interactions with the environment.
  fit <- labraduck(Y = Y, upsilon = taxa_covariance$upsilon, Xi = taxa_covariance$Xi, F = F, G = G, W = W, M0 = M0, C0 = C0,
                   observations = days, gamma_scale = (var_scale * 2/3), W_scale = (var_scale * 1/3), apply_smoother = MAP,
                   n_samples = n_samples, ret_mean = MAP)

  # pack up results
  fit_obj <- list(Y = Y, alr_ys = alr_ys, alr_ref = alr_ref, fit = fit)

  # save results
  if(MAP) {
    save_dir <- check_output_dir(c("output","model_fits",paste0(tax_level,"_MAP")))
  } else {
    save_dir <- check_output_dir(c("output","model_fits",tax_level))
  }
  saveRDS(fit_obj, file.path(save_dir,paste0(host,"_labraduckfit.rds")))
}

#' Perform k-fold cross-validation for a given host and kernel configuration choice
#' 
#' @param data a phyloseq object
#' @param host host short name (e.g. ACA)
#' @param sample_covariance composite kernel function
#' @param rho scalar bandwidth parameter
#' @param tax_level taxonomic level at which to agglomerate data
#' @param holdout_proportion proportion of host's sample to use as a test set
#' @param iterations number of CV iterations to perform; this should be <= n_samples for this host
#' @details average error (currently RMSE of log counts)
#' @return NULL
#' @import phyloseq
#' @export
#' @examples
#' tax_level <- "ASV"
#' data <- load_data(tax_level = tax_level)
#' sample_covariance <- get_Gamma(kernel_scale = 2, proportions = c(1, 0, 0), min_correlation = 0.1, days_to_baseline = 90)
#' taxa_covariance <- get_Xi(phyloseq::ntaxa(data), total_variance = 1)
#' perform_cv(data, host = "GAB", taxa_covariance = taxa_covariance, sample_covariance = sample_covariance, tax_level = tax_level)
perform_cv <- function(data, host, taxa_covariance, sample_covariance, tax_level = "ASV", holdout_proportion = 0.2) {
  params <- formalize_parameters(data)
  return(fit_GP(data, host, taxa_covariance, sample_covariance, rho, tax_level, alr_ref = params$alr_ref, n_samples = 100,
                             MAP = FALSE, holdout_proportion = holdout_proportion))
}








