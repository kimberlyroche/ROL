# quick script to calculate autocorrelation

library(ROL)

tax_level <- "ASV"
data <- load_data(tax_level="ASV", host_sample_min=75, count_threshold=5, sample_threshold=0.2)
lagged_ac <- calc_autocorrelation(data, lag_units="months", lag_max=60, use_lr="ilr", alr_ref=NULL, resample=TRUE, resample_rate=0.2)
plot_mean_autocorrelation(lagged_ac, show_plot=FALSE)
