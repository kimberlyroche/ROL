library(ROL)
library(phyloseq)
library(filelock)

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1) {
  stop("Missing argument for host!\n")
}
host <- args[1]

record_error <- function(results, label) {
  out_file <- file.path("output","cv_results.txt")
  lock_file <- paste0(out_file, ".lck")
  lck <- lock(lock_file, timeout = Inf)
  write(paste0(results$fold_error,"\tfold_error\t",host,"\t",label), file = out_file, append = TRUE)
  write(paste0(results$log_rmse,"\tlog rmse\t",host,"\t",label), file = out_file, append = TRUE)
  unlock(lck)
}

tax_level <- "ASV"
data <- load_data(tax_level="ASV", host_sample_min=75, count_threshold=5, sample_threshold=0.2)

# evaluate 3 models:
# (1) SE kernel only, with large bandwidth (autocorrelation decays away at 90 days)
# (2) SE kernel only, with small bandwidth ("flexible"; autocorrelation decays away at 30 days)
# (3) SE kernel + seasonal components (incl. diet and climate)

rho_30 <- calc_se_decay(min_correlation = 0.1, days_to_baseline = 30)
rho_90 <- calc_se_decay(min_correlation = 0.1, days_to_baseline = 90)

Gamma.1 <- get_Gamma(kernel_scale = 2, proportions = c(1, 0, 0), rho = rho_90)
Gamma.2 <- get_Gamma(kernel_scale = 2, proportions = c(1, 0, 0), rho = rho_30)
Gamma.3 <- get_Gamma(kernel_scale = 2, proportions = c(0.5, 0.25, 0.25), rho = rho_90)

cat("Fitting model 1...\n")
e1 <- perform_cv(data, host = host, composite_kernel = Gamma.1, tax_level = tax_level, holdout_proportion = 0.2)
record_error(e1, "SE_wide")

cat("Fitting model 2...\n")
e2 <- perform_cv(data, host = host, composite_kernel = Gamma.2, tax_level = tax_level, holdout_proportion = 0.2)
record_error(e2, "SE_narrow")

cat("Fitting model 3...\n")
e3 <- perform_cv(data, host = host, composite_kernel = Gamma.3, tax_level = tax_level, holdout_proportion = 0.2)
record_error(e3, "seasonal")

# calculate autocorrelation; previous
# lagged_ac <- calc_autocorrelation(data, lag_units="months", lag_max=60, use_lr="ilr", alr_ref=NULL, resample=TRUE, resample_rate=0.2)
# plot_mean_autocorrelation(lagged_ac, show_plot=FALSE)
