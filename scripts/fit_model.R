library(ROL)
library(phyloseq)
library(fido)
library(optparse)

option_list = list(
  make_option(c("--host"), type = "character", default = NULL,
              help = "host short name", metavar = "character"),
  make_option(c("--tax_level"), type = "character", default = NULL,
              help = "taxonomic level", metavar = "character"),
  make_option(c("--MAP"), type = "logical", default = FALSE,
              help = "use MAP fit (point estimate)", metavar = "logical"),
  make_option(c("--DLM"), type = "logical", default = FALSE,
              help = "use DLM", metavar = "logical"),
  make_option(c("--cov"), type = "logical", default = TRUE,
              help = "use covariates", metavar = "logical")
  # make_option(c("--scramble"), type = "logical", default = FALSE,
  #             help = "scramble taxa and time points to simulate null", metavar = "logical")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$host) | is.null(opt$tax_level)) {
  stop("Missing host and taxonomic level!")
}

data <- load_data(tax_level = opt$tax_level, host_sample_min = 75, count_threshold = 5, sample_threshold = 0.2)
params <- formalize_parameters(data)

if(opt$DLM) {
  taxa_covariance <- get_Xi(phyloseq::ntaxa(data), total_variance = 1)
  fit_DLM(data, host = opt$host, taxa_covariance, var_scale = 1, tax_level = opt$tax_level, alr_ref = params$alr_ref, MAP = opt$MAP, use_covariates = opt$cov)
} else {
  # use GP
  sample_covariance <- get_Gamma(kernel_scale = 2, proportions = c(1/3, 1/3, 1/3), min_correlation = 0.1, days_to_baseline = 90)
  taxa_covariance <- get_Xi(phyloseq::ntaxa(data), total_variance = 1)
  fit_GP(data, host = opt$host, taxa_covariance = taxa_covariance, sample_covariance = sample_covariance, tax_level = opt$tax_level, alr_ref = params$alr_ref, MAP = opt$MAP)
}
