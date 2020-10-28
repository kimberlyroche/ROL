# Compare model fits with and without covariates in terms of (1) "rug" plots and (2) universality of associations
# Assumes no-covariate MAP fits are in the directory: output/model_fits_nocovariates/ASV_MAP
# Assumes covariate-inclusive MAP fits are in the directory: output/model_fits/ASV_MAP

library(fido)
library(ROL)
library(phyloseq)
library(dplyr)
library(driver)
library(ggplot2)

# Read covariate-free MAP model fits
get_Sigmas_DLM <- function(model_dir) {
  pattern_str <- "*_labraduckfit.rds"
  regexpr_str <- "_labraduckfit.rds"
  level_dir <- paste0("output/model_fits/ASV_MAP/", model_dir)
  model_list <- list.files(path = level_dir, pattern = pattern_str, full.names = TRUE, recursive = FALSE)
  hosts <- as.vector(sapply(model_list, function(x) { idx <- regexpr(regexpr_str, x); return(substr(x, idx-3, idx-1)) } ))
  Sigmas <- list()
  for(i in 1:length(hosts)) {
    host <- hosts[i]
    cat("Parsing host",i,"/",length(hosts)," (",host,")\n")
    fit <- readRDS(model_list[i])$fit
    fit <- to_clr(fit)
    Sigmas[[host]] <- fit$Sigma[,,1]
  }
  return(Sigmas)
}

Sigmas_cov <- get_Sigmas_DLM(model_dir = "covariates")
Sigmas_nocov <- get_Sigmas_DLM(model_dir = "no_covariates")
# Sigmas_noise <- get_Sigmas_DLM(model_dir = "noise_covariates")

coordinate_number <- dim(Sigmas_cov[[1]])[1]
# label_pairs <- matrix(NA, coordinate_number, coordinate_number)
# for(i in 1:coordinate_number) {
#   for(j in 1:coordinate_number) {
#     if(i < j) {
#       label_pairs[i,j] <- paste0(i,"_",j)
#     }
#   }
# }
# label_pairs <- label_pairs[upper.tri(label_pairs, diag=F)]

get_association_matrix <- function(Sigmas) {
  hosts <- names(Sigmas)
  H <- length(hosts)
  associations <- matrix(NA, H, (coordinate_number^2)/2 - coordinate_number/2)
  for(m in 1:H) {
    host <- hosts[m]
    Sigma <- Sigmas[[m]]
    rhos <- c()
    for(j in 2:coordinate_number) {
      for(i in 1:(j-1)) {
        var_i <- Sigma[i,i]
        var_j <- Sigma[j,j]
        var_i_minus_j <- var_i + var_j - 2*Sigma[i,j]
        rho_ij <- 1 - var_i_minus_j / (var_i + var_j)
        rhos <- c(rhos, rho_ij)
      }
    }
    associations[m,] <- rhos
  }
  associations
}

assoc_cov <- get_association_matrix(Sigmas_cov)
assoc_nocov <- get_association_matrix(Sigmas_nocov)
# assoc_noise <- get_association_matrix(Sigmas_noise)

# Calculate a canonical column and row ordering
d.r <- dist(assoc_cov)
o.r <- hclust(d.r)$order
d.c <- dist(t(assoc_cov))
o.c <- hclust(d.c)$order

assoc_cov_clustered <- assoc_cov[o.r,o.c]
assoc_nocov_clustered <- assoc_nocov[o.r,o.c]
# assoc_noise_clustered <- assoc_noise[o.r,o.c]
# label_pairs_clustered <- label_pairs[o.c]

render_heatmap <- function(associations, save_file) {
  df <- gather_array(associations, "proportionality", "host", "pair")
  p <- ggplot(df, aes(pair, host)) +
    geom_tile(aes(fill = proportionality)) +
    scale_fill_gradient2(low = "darkblue", high = "darkred")
  ggsave(save_file, p, units="in", dpi=150, height=5, width=15)
}

render_heatmap(assoc_cov_clustered, "rug_cov.png")
render_heatmap(assoc_nocov_clustered, "rug_nocov.png")
# render_heatmap(assoc_noise_clustered, "rug_noise.png")
