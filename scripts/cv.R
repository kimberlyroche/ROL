library(ROL)
library(phyloseq)
library(filelock)

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1) {
  stop("Missing argument for host!\n")
}
host <- args[1]

record_error <- function(results, days_to_baseline, incl_covariates, Gamma_scale, Xi_scale) {
  out_file <- file.path("output","cv_results.txt")
  lock_file <- paste0(out_file, ".lck")
  lck <- lock(lock_file, timeout = Inf)
  write(paste0(round(results$fold_error,5),"\tfold_error\t",host,"\t",days_to_baseline,"\t",
               incl_covariates,"\t",Gamma_scale,"\t",Xi_scale), file = out_file, append = TRUE)
  write(paste0(round(results$log_rmse,5),"\tlog rmse\t",host,"\t",days_to_baseline,"\t",
               incl_covariates,"\t",Gamma_scale,"\t",Xi_scale), file = out_file, append = TRUE)
  unlock(lck)
}

record_failed_fit <- function(out_str) {
  out_file <- file.path("output","cv_failed_fits.txt")
  lock_file <- paste0(out_file, ".lck")
  lck <- lock(lock_file, timeout = Inf)
  write(out_str, file = out_file, append = TRUE)
  unlock(lck)
}

tax_level <- "ASV"
data <- load_data(tax_level="ASV", host_sample_min=75, count_threshold=5, sample_threshold=0.2)
params <- formalize_parameters(data)

# cycle through evaluation criteria:
# (1) Gamma: large x small bandwidth
# (2)        include x exclude covariates
# (3)        large x small total variance
# (4)    Xi: large x small total variance

for(days_to_baseline in c(90, 30)) {
  for(incl_covariates in c(TRUE, FALSE)) {
    for(Gamma_scale in c(2, 1)) {
      for(Xi_scale in c(2, 1)) {
        cat(paste0("Evaluating ",host," on DTB=",days_to_baseline,", COV=",incl_covariates,", GS=",Gamma_scale,", XS=",Xi_scale,"\n"))
        if(incl_covariates) {
          kernel_proportions <- c(0.5, 0.25, 0.25)
        } else {
          kernel_proportions <- c(1, 0, 0)
        }
        sample_covariance <- get_Gamma(kernel_scale = Gamma_scale, proportions = kernel_proportions,
                               min_correlation = 0.1, days_to_baseline = days_to_baseline)
        taxa_covariance <- get_Xi(phyloseq::ntaxa(data), total_variance = Xi_scale)
        # this can bomb if model fit fails to converge during optimization; for now, just pass over this case
        # we need to think about whether this is systematically affecting error in these cases
        err = tryCatch({
          perform_cv(data, host, taxa_covariance, sample_covariance, tax_level, holdout_proportion = 0.2)
        }, error = function(e) {
          # write to file
          out_str <- paste0(host,"\t",days_to_baseline,"\t",incl_covariates,"\t",Gamma_scale,"\t",Xi_scale)
          record_failed_fit(out_str)
        }, finally = {})
        if(!is.null(err)) {
          record_error(err, days_to_baseline, incl_covariates, Gamma_scale, Xi_scale)
        }
      }
    }
  }
}
