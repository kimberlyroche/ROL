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
#' @param MAP use MAP estimate model output instead of full posterior output
#' @return NULL
#' @export
#' @examples
#' model_list <- get_fitted_model_list(tax_level="ASV", MAP=FALSE)
get_fitted_model_list <- function(tax_level="ASV", MAP=FALSE) {
  # all fitted models have this suffix
  pattern_str <- "*_bassetfit.rds"
  regexpr_str <- "_bassetfit.rds"
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
#' @param log_var_scale scale of the log variance
#' @return list containing inverse Wishart parameters degrees of freedom and scale matrix
#' @export
#' @examples
#' params <- default_ALR_prior(D=100, log_var_scale=1)
#' Sigma <- matrixsampling::rinvwishart(1, params$upsilon, params$Xi)[,,1]
default_ALR_prior <- function(D, log_var_scale=1) {
  upsilon <- D-1+10 # specify low certainty/concentration
  GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference
  Xi <- GG%*%(diag(D)*log_var_scale)%*%t(GG) # take diag as covariance over log abundances
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

#' Fit a Gaussian process to a single host series using basset
#' 
#' @param data a phyloseq object
#' @param host host short name (e.g. ACA)
#' @param tax_level taxonomic level at which to agglomerate data
#' @param SE_days_to_baseline days at which squared exponential kernel decays to baseline correlation
#' @param date_lower_limit minimum date to consider (string format: YYYY-MM-DD)
#' @param date_upper_limit maximum date to consider (string format: YYYY-MM-DD)
#' @param alr_ref index of reference ALR coordinate
#' @param n_samples number of posterior samples to draw
#' @param MAP compute MAP estimate only (as single posterior sample)
#' @details Fitted model and metadata saved to designated model output directory.
#' @return NULL
#' @import phyloseq
#' @import stray
#' @export
#' @examples
#' data <- load_data(tax_level="ASV")
#' fit_GP(data, host="ACA", tax_level="ASV", SE_sigma=1, PER_sigma=1, SE_days_to_baseline=90, MAP=FALSE)
fit_GP <- function(data, host, tax_level="ASV", SE_days_to_baseline=90, date_lower_limit=NULL, date_upper_limit=NULL,
                   alr_ref=NULL, n_samples=100, MAP=FALSE, ...) {
  if(MAP) {
    cat(paste0("Fitting stray::basset model (MAP) to host ",host,"\n"))
  } else {
    cat(paste0("Fitting stray::basset model to host ",host,"\n"))
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
  host.num <- as.character(data.name_mapping[data.name_mapping$sname == host,]$host_subject_id2)
  data.diet <- subset_samples(data.diet, host == host.num)
  metadata.diet <- sample_data(data.diet)

  baseline_date <- host_metadata$collection_date[1]
  observations <- sapply(host_metadata$collection_date, function(x) round(difftime(x, baseline_date, units="days"))) + 1
  
  # pull out the count table
  Y <- otu_table(host_data)@.Data

  X <- build_design_matrix(host_metadata, metadata.diet)
  
  # if a specified data lower and upper limit have been specified, chop the observed series down
  # to this range
  lower_idx <- NULL
  upper_idx <- NULL
  if(!is.null(date_lower_limit) & !is.null(date_upper_limit)) {
    lower_idx <- which(names(observations) >= date_lower_limit)
    upper_idx <- which(names(observations) <= date_upper_limit)
  }
  if(length(lower_idx) > 0 & length(upper_idx) > 0) {
    unique_lower_idx <- min(lower_idx)
    unique_upper_idx <- max(upper_idx)
    # subset dates
    Y_pre <- Y
    Y <- t(Y_pre[unique_lower_idx:unique_upper_idx,])
    rm(Y_pre)
    # set first observation to t=1
    observations_pre <- observations
    observations <- matrix(observations_pre[unique_lower_idx:unique_upper_idx], nrow=1) - observations_pre[unique_lower_idx] + 1
    rm(observations_pre)
  } else {
    # just clean up
    Y <- t(Y)
    dim(observations) <- c(1, length(observations))
  }
  
  # get dimensions
  D <- nrow(Y)
  N <- ncol(Y)
  
  # stray uses the D^th element as the ALR reference by default
  # if we'd like to use a different reference, do some row shuffling in Y to put the reference at the end
  if(!is.null(alr_ref)) {
    Y <- Y[c(setdiff(1:D,alr_ref),alr_ref),]
  }
  
  # strip off sequence variant labels
  colnames(Y) <- NULL
  rownames(Y) <- NULL
  
  # back-calculate the squared exponential bandwidth parameter by finding a bandwidth that gives
  # a desired minimum correlation at the number of days specified by SE_days_to_baseline  
  dc <- 0.1 # desired minimum correlation
  rho <- sqrt(-SE_days_to_baseline^2/(2*log(dc))) # back calculate the decay
  
  # define the composite kernel over samples
  Gamma <- function(X) {
    jitter <- 1e-08
    # Gamma_scale <- params$total_variance
    Gamma_scale <- 2
    base_sigma_sq <- Gamma_scale * 0.5
    PER_sigma_sq <- Gamma_scale * 0.5
    SE(X[1,,drop=F], sigma = sqrt(base_sigma_sq), rho = rho, jitter = jitter) +
      PER(X[1,,drop=F], sigma = sqrt(PER_sigma_sq/2), rho = 1, period = 365, jitter = jitter) +
      SE(X[2:6,,drop=F], sigma = sqrt(PER_sigma_sq/2), rho = 1, jitter = jitter) # a few diet PCs
  }
  
  # use a default ALR prior
  prior_params <- default_ALR_prior(D)
  
  # define the prior over baselines
  alr_ys <- driver::alr((t(Y) + 0.5))
  alr_means <- colMeans(alr_ys)
  Theta <- function(X) matrix(alr_means, D-1, ncol(X))
  
  if(MAP) {
    n_samples <- 0
  }
  fit <- stray::basset(Y, X, prior_params$upsilon, Theta, Gamma, prior_params$Xi,
                       n_samples = n_samples, ret_mean = MAP, 
                       b2 = 0.98, step_size = 0.004, eps_f = 1e-11, eps_g = 1e-05,
                       max_iter = 10000L, optim_method = "adam")

  if(MAP) {
    # fill out dimensions; some function (including stray predictive functions) expect arrays
    # not matrices here
    fit <- fix_MAP_dims(fit)
  }

  # pack up results
  fit_obj <- list(Y = Y, alr_ys = alr_ys, alr_ref = alr_ref, X = X, fit = fit)
  
  # a hack: for later prediction, these will need to be available in the workspace for Gamma
  fit_obj$kernelparams$SE_rho <- rho

  # save results
  if(MAP) {
    save_dir <- check_output_dir(c("output","model_fits",paste0(tax_level,"_MAP")))
  } else {
    save_dir <- check_output_dir(c("output","model_fits",tax_level))
  }
  saveRDS(fit_obj, file.path(save_dir,paste0(host,"_bassetfit.rds")))
}
