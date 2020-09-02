# JLM visualizations

library(ROL)
library(ggplot2)

# (1) all utilized samples across time

data <- load_data(tax_level = "ASV")
md <- sample_data(data)

plot_df <- data.frame(sample = md$collection_date, host = md$sname)
head(plot_df)

ggplot(plot_df, aes(x = sample, y = host)) +
  geom_point()

# (2) binnned correlation between hosts

Sigmas <- load_MAP_estimates(tax_level = "ASV", logratio = "clr")
interactions <- get_pairwise_correlations(tax_level = "ASV", logratio = "clr", Sigmas = Sigmas)
str(interactions)

DUI_idx <- which(names(Sigmas) == "DUI")
reference_correlations <- abs(interactions$interactions[DUI_idx,])
n_hosts <- nrow(interactions$interactions)
max_sd <- sd(c(rep(-1, ceiling(n_hosts/2)), rep(1, floor(n_hosts)/2)))

# plot(density(reference_correlations))
# these boundaries for "low", "medium", "high" seem reasonable
boundaries <- c(0.3, 0.5, 1.0)
boundaries <- c(0, boundaries)
binned_list <- list()
for(b in 2:length(boundaries)) {
  idx <- which(reference_correlations > boundaries[b-1] & reference_correlations <= boundaries[b])
  subset_interactions <- interactions$interactions[,idx]
  sds <- apply(subset_interactions, 2, sd)
  binned_list[[b-1]] <- sds
}

plot_df <- data.frame(obs_sd = binned_list[[1]], quantile = "low")
plot_df <- rbind(plot_df, data.frame(obs_sd = binned_list[[2]], quantile = "medium"))
plot_df <- rbind(plot_df, data.frame(obs_sd = binned_list[[3]], quantile = "high"))
plot_df$quantile <- as.factor(plot_df$quantile)
# fix levels
levels(plot_df$quantile) <- levels(plot_df$quantile)[c(2,3,1)]

ggplot(plot_df, aes(x = obs_sd, color = quantile)) +
  geom_density(size = 2) +
  xlab("standard deviation across hosts") +
  ylab("density")

