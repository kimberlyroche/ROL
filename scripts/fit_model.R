library(ROL)
library(phyloseq)
library(stray)

args = commandArgs(trailingOnly=TRUE)

host <- args[1]
tax_level <- args[2]
MAP <- as.logical(args[3])
scramble <- as.logical(args[4])

data <- load_data(tax_level=tax_level, host_sample_min=75, count_threshold=5, sample_threshold=0.2)
params <- formalize_parameters(data)

# no covariate model
# sample_covariance <- get_Gamma(kernel_scale = 2, proportions = c(1, 0, 0), min_correlation = 0.1, days_to_baseline = 90)
# equal contribution of autocorrelation, season, and diet and rainfall PCs
sample_covariance <- get_Gamma(kernel_scale = 2, proportions = c(1/3, 1/3, 1/3), min_correlation = 0.1, days_to_baseline = 90)
taxa_covariance <- get_Xi(phyloseq::ntaxa(data), total_variance = 1)

fit_GP(data, host = host, taxa_covariance = taxa_covariance, sample_covariance = sample_covariance, tax_level = tax_level, alr_ref = params$alr_ref, MAP = MAP, scramble = scramble)
