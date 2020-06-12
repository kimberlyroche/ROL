library(ROL)
library(phyloseq)

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1) {
  stop("Missing argument for host!\n")
}
host <- args[1]

tax_level <- "ASV"
data <- load_data(tax_level="ASV", host_sample_min=75, count_threshold=5, sample_threshold=0.2)

Gamma.se_only <- get_Gamma(kernel_scale = 2, proportions = c(1, 0, 0), days_to_baseline = 90)
Gamma.seasonal <- get_Gamma(kernel_scale = 2, proportions = c(0.5, 0.25, 0.25), days_to_baseline = 90)
out_file <- file.path("output","cv_results.txt")

e1 <- perform_cv(data, host = host, kernels = Gamma.se_only, tax_level = tax_level, holdout_proportion = 0.2)
write(paste0(e1,"\tlog rmse\t",host,"\tSE_only"), file = out_file, append = TRUE)
e2 <- perform_cv(data, host = host, kernels = Gamma.seasonal, tax_level = tax_level, holdout_proportion = 0.2)
write(paste0(e2,"\tlog rmse\t",host,"\tseasonal"), file = out_file, append = TRUE)

# calculate autocorrelation
# lagged_ac <- calc_autocorrelation(data, lag_units="months", lag_max=60, use_lr="ilr", alr_ref=NULL, resample=TRUE, resample_rate=0.2)
# plot_mean_autocorrelation(lagged_ac, show_plot=FALSE)
