# Fit GP to each host's PC1 in ABRP data set
# Find a per-host fit that reduces lag2 AC < 0.2
# Use the residual of this model fit as the "detrended" PC1 and append it to metadata

library(kernlab)
library(phyloseq)
library(ggplot2)
library(gridExtra)

# ----------------------------------------------------------------------------------------------
#   PARSE DATA
# ----------------------------------------------------------------------------------------------

# load the version of the data set with host short names
data <- readRDS("input/ps0.RDS")
full_metadata <- sample_data(data)
host_list <- unique(full_metadata$sname)

# load the version of the data set with PCs in the metadata
data <- readRDS("input/ps.RDS")
full_metadata <- sample_data(data)

# load the host label mapping: "Baboon_###" to "XXX"
mapping <- read.table("input/host_subject_id_to_sname_key.csv", sep = ",", header = TRUE)

# tack on an extra column to the metadata
full_metadata$modified_PC1 <- numeric(phyloseq::nsamples(data))

# ----------------------------------------------------------------------------------------------
#   VARIOUS FUNCTIONS
# ----------------------------------------------------------------------------------------------

# function to pad observations with NAs
# note: in the rare case of duplicate samples for a given host for a given day
#       (replicates?), I'm placing those samples next to each other in the observation
#       vector, e.g.
#            1  NA  NA  NA  5  5  NA  NA  8  NA  NA  10  11  NA  12  ...
#                           ^  ^
#                        duplicates
pad_observations <- function(x, y) {
  n_duplicates <- length(x) - length(unique(x))
  x_padded <- rep(NA, max(x) + n_duplicates)
  offset <- 0
  for(xx in x) {
    if(!is.na(x_padded[xx + offset])) {
      offset <- offset + 1
    }
    x_padded[xx + offset] <- xx
  }
  y_padded <- x_padded
  y_padded[!is.na(y_padded)] <- y
  return(list(x = x_padded, y = y_padded))
}

# remove autocorrelation by fitting either a Gaussian process (via gausspr) or
# an AR(1) model (via arima) and saving the residual
# note: because of the irregularities in sampling frequency, the model tends to fit
#       individual host series with varying amounts of flexibility; in general
#       hard-coding the parameters (e.g. length scale) across hosts didn't seem to
#       improve the average fit, either by cross-validation or by reduction in the AcF
#       so for now I'm letting the model auto-fit parameters
detrend <- function(subset_data, which_method = "GP", visualize = FALSE) {
  metadata <- sample_data(subset_data)
  x <- metadata$collection_date
  y <- metadata$PC1
  PC1_baseline <- mean(y)
  y <- scale(y, center = TRUE, scale = FALSE)
  baseline_date <- min(x)
  x <- sapply(x, function(z) difftime(z, baseline_date)) + 1
  duplicates <- as.numeric(names(which(table(x) > 1)))
  
  # plot residual autocorrelation
  padded_series <- pad_observations(x, y)
  x <- padded_series$x
  y <- padded_series$y
  
  if(which_method == "GP") {
    # fit GP
    fit <- gausspr(x, y, scaled = TRUE)
    # the models below seem to overfit, as does the AR(1) model generally
    # both low var and high sigma give the model a lot of flexibility
    # fit <- gausspr(x, y, scaled = TRUE, var = 0.002)
    # fit <- gausspr(x, y, scaled = TRUE, kpar = list(sigma = 10), var = 0.01)
    mean_prediction <- predict(fit, x)
    residual <- y - mean_prediction + PC1_baseline
  } else {
    # AR(1)
    # this seems to work OK but there are occasional convergence issues;
    # exhaustively trying all available optimization methods seems to be a reasonable
    # workaround, which is what I'm doing below
    fit_methods <- c("L-BFGS-B", "BFGS", "Nelder-Mead", "CG", "SANN", "Brent")
    for(fit_method in fit_methods) {
      fit <- tryCatch({ arima(y, order = c(1, 0, 0), optim.method = fit_method) },
                      warning = function(w) { },
                      error = function(e) { NULL },
                      finally = {})
      if(!is.null(fit)) {
        break
      }
    }
    if(is.null(fit)) {
      return(list(residuals = NULL, lag_acf = NULL))
    }
    mean_prediction <- y - fit$residuals
    residual <- fit$residuals + PC1 + baseline
  }
  
  if(visualize) {
    # visualize (1) model fit and (2) AC before (3) AC after (residual)
    par(mfrow = c(2,2))
    plot(x, y, type ="p")
    lines(x[!is.na(x)], mean_prediction[!is.na(x)], type = "l", col = "red")
    plot(x, residual, type ="p")
    acf(y, na.action = na.pass)
    acf(residual, na.action = na.pass)
  }

  acf_data <- acf(residual, na.action = na.pass, plot = FALSE)
  return(list(residual = residual, lag_acf = acf_data$acf[,1,1], duplicates = duplicates))
}

# get the observed autocorrelation associated with the shortest available lag
# (e.g. 2 days, 3 days, ... 10 days, etc.)
get_min_lag <- function(lags) {
  available_lags <- which(!is.na(lags))
  available_lags <- available_lags[2:length(available_lags)]
  report_lag <- min(available_lags)
  return(lags[report_lag])
}

# ----------------------------------------------------------------------------------------------
#   FIT THE GP AND EXTRACT THE RESIDUALS
# ----------------------------------------------------------------------------------------------

original_acf <- c()
modified_acf <- c()
host_sample_size <- c()
excluded_from_fit <- c()
visualize_flag <- TRUE # we'll visualize the first individual's fit, then skip the rest
for(host in host_list) {
  cat("Fitting host:",host,"\n")
  # translate to a label of the form "Baboon_###"
  host_id <- unique(mapping[mapping$sname == host,]$host_subject_id2)
  subset_data <- subset_samples(data, host == host_id)
  host_sample_n <- phyloseq::nsamples(subset_data)
  # this commented out code checks to see if the collection times are in order; we're
  # assuming they are when we fill in the modified PC values
  # they do appear to be in order for all animals
  #   host_dates <- full_metadata[full_metadata$host == host_id,]$collection_date
  #   date_order_compare <- order(host_dates) != order(sort(host_dates))
  #   if(sum(date_order_compare) > 0) {
  #     cat("\tReplacement indices are not in order!\n")
  #   }
  # get the indices (in the metagenomics data set) of the PC values to replace
  # for this host
  replace_idx <- which(full_metadata$host == host_id)
  if(host_sample_n > 10) {
    result.GP <- detrend(subset_data, which_method = "GP", visualize = visualize_flag)
    visualize_flag <- FALSE
    full_metadata$modified_PC1[replace_idx] <- result.GP$residual[!is.na(result.GP$residual)]
    original_acf <- c(original_acf,
                      get_min_lag(acf(full_metadata$PC1[replace_idx], na.action = na.pass, plot = FALSE)$acf[,1,1]))
    modified_acf <- c(modified_acf,
                      get_min_lag(result.GP$lag_acf))
    host_sample_size <- c(host_sample_size, host_sample_n)
    excluded_from_fit <- c(excluded_from_fit, rep(FALSE, host_sample_n))
  } else {
    # if this host has fewer than 10 samples, skip the model fitting;
    # I'm reasonably certain AC won't be a problem for these individuals
    full_metadata$modified_PC1[replace_idx] <- full_metadata$PC1[replace_idx]
    excluded_from_fit <- c(excluded_from_fit, rep(TRUE, host_sample_n))
  }
}

# ----------------------------------------------------------------------------------------------
#   VISUALIZATION
# ----------------------------------------------------------------------------------------------

# (1) Show the overall reduction in autocorrelation associated with the shortest lag for each host
plot_df <- data.frame(x = original_acf, label = "ORIGINAL DATA")
plot_df <- rbind(plot_df, data.frame(x = modified_acf, label = "MODIFIED DATA"))
ggplot(plot_df, aes(x = x)) +
  geom_histogram(bins = 30) +
  facet_grid(rows = vars(label)) +
  xlab("ACF") +
  xlim(-1, 1)

# (2) Show all values associated with the original and modified PC1; there's still an obvious
#     seasonal signal, at least at *scale*
selected_host <- "Baboon_4"
plot_df <- data.frame(x = full_metadata$collection_date,
                      y = full_metadata$PC1,
                      label = "ORIGINAL DATA",
                      label2 = full_metadata$host == selected_host,
                      host_type = as.factor(excluded_from_fit))
plot_df <- rbind(plot_df, data.frame(x = full_metadata$collection_date,
                                     y = full_metadata$modified_PC1,
                                     label = "MODIFIED DATA",
                                     label2 = full_metadata$host == selected_host,
                                     host_type = as.factor(excluded_from_fit)))
p1 <- ggplot() +
  geom_point(data = plot_df[plot_df$label == "ORIGINAL DATA" & plot_df$label2 == FALSE,], aes(x = x, y = y), color = "#AAAAAA") +
  geom_point(data = plot_df[plot_df$label == "ORIGINAL DATA" & plot_df$label2 == TRUE,], aes(x = x, y = y), color = "blue") +
  xlab("time index") +
  ylab("PC1")
p2 <- ggplot() +
  geom_point(data = plot_df[plot_df$label == "MODIFIED DATA" & plot_df$label2 == FALSE,], aes(x = x, y = y), color = "#AAAAAA") +
  geom_point(data = plot_df[plot_df$label == "MODIFIED DATA" & plot_df$label2 == TRUE,], aes(x = x, y = y), color = "blue") +
  xlab("time index") +
  ylab("PC1")
grid.arrange(p1, p2, nrow = 2)

# ----------------------------------------------------------------------------------------------
#   SAVE OUTPUT
# ----------------------------------------------------------------------------------------------

sample_data(data) <- full_metadata
saveRDS(data, file = "input/ps_modifiedPC1.rds")


