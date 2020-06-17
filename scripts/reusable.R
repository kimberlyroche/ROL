library(ROL)

calculate autocorrelation; previous
lagged_ac <- calc_autocorrelation(data, lag_units="months", lag_max=60, use_lr="ilr", alr_ref=NULL, resample=TRUE, resample_rate=0.2)
plot_mean_autocorrelation(lagged_ac, show_plot=FALSE)
