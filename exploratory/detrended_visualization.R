library(phyloseq)
library(ROL)
library(ggplot2)

# to do: (1) run this for many different individuals
#        (2) incorporate their lagged AC
#        (3) ribbon plot for Lambda, Lambda_detrended

tax_level <- "ASV"
data <- load_data(tax_level = tax_level, host_sample_min = 75, count_threshold = 5, sample_threshold = 0.2)
metadata <- sample_data(data)
cat(paste0("There are ",length(unique(metadata$sname))," unique hosts in this data set and a total of ",phyloseq::nsamples(data)," samples!\n"))
params <- formalize_parameters(data)

for(host in unique(metadata$sname)) {
  cat("Fitting",host,"\n")
  fit_GP(data, host = host, tax_level = tax_level, SE_days_to_baseline = 90, alr_ref = params$alr_ref, MAP = TRUE)
}

use_Eta <- FALSE
use_detrended <- FALSE

df <- data.frame(x = c(), y = c(), host = c())
for(host in unique(metadata$sname)) {
  fit <- readRDS(paste0("output/model_fits/ASV_MAP/",host,"_bassetfit.rds"))

  # calculate autocorrelation for Lambda
  observations <- fit$X[1,]
  Lambda <- fit$fit$Lambda[,,1] # this is D-1 taxa x N samples
  Eta <- fit$fit$Eta[,,1]
  Theta <- fit$fit$Theta(fit$fit$X)
  Gamma <- fit$fit$Gamma(fit$fit$X)
  Gamma_sqrt <- chol(Gamma)
  
  Lambda_detrended <- Theta + (Lambda - Theta)%*%solve(Gamma_sqrt)

  lags <- list()
  for(i in 1:(length(observations)-1)) {
    for(j in (i+1):length(observations)) {
      diff_week <- round(abs(observations[i] - observations[j])/7)
      if(is.na(diff_week)) {
        cat(i,",",j,"\n")
      }
      lag_str <- as.character(diff_week)
      if(use_Eta) {
        if(use_detrended) {
          Eta_detrended <- Lambda_detrended + (Eta - Lambda)
          ij_correlation <- cor(Eta_detrended[,i], Eta_detrended[,j])
        } else {
          ij_correlation <- cor(fit$fit$Eta[,i,1], fit$fit$Eta[,j,1])
        }
      } else {
        if(use_detrended) {
          ij_correlation <- cor(Lambda_detrended[,i], Lambda_detrended[,j])
        } else {
          ij_correlation <- cor(Lambda[,i], Lambda[,j])
        }
      }
      if(lag_str %in% names(lags)) {
        lags[[lag_str]] <- c(lags[[lag_str]], ij_correlation)
      } else {
        lags[[lag_str]] <- c(ij_correlation)
      }
    }
  }
  df <- rbind(df, data.frame(x = as.numeric(names(lags)), y = sapply(lags, function(lag) mean(lag)), host = host))
}

p <- ggplot(df[df$x < 104,]) +
  geom_smooth(aes(x = x, y = y)) +
  geom_point(aes(x = x, y = y)) +
  xlab("lag (weeks)") +
  ylab("ACF")
show(p)
if(use_Eta) {
  save_file <- "Eta"
} else {
  save_file <- "Lambda"
}
if(use_detrended) {
  save_file <- paste0(save_file, "_detrended")
}
save_file <- paste0(save_file, ".png")
ggsave(paste0("C:/Users/kim/Desktop/",save_file), p, dpi = 100, units = "in", height = 6, width = 10)






