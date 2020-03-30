library(ggplot2)
library(gridExtra)
library(driver)

# get distributions of distances between: (1) "baselines": empirical average ALR abundances
#                                         (2) "dynamics": microbial covariance matrices
model_list <- get_fitted_model_list(tax_level="ASV", MAP=TRUE)
baselines <- matrix(NA, model_list$D-1, length(model_list$model_list))
for(i in 1:length(model_list$model_list)) {
  data <- readRDS(model_list$model_list[i])
  baselines[,i] <- colMeans(data$alr_ys)
}
baseline_dist <- as.matrix(dist(t(baselines)))
dynamics_dist <- readRDS("output/Sigma_distance_ASV_MAP.rds")$distance_mat

# order of hosts in baseline and dynamics distance matrices is demonstrably
# identical as in model_list$hosts

baseline_distro <- baseline_dist[upper.tri(baseline_dist, diag=F)]
dynamics_distro <- dynamics_dist[upper.tri(dynamics_dist, diag=F)]

# visualize distributions over these distances
if(FALSE) {
  df.baseline <- data.frame(distance=c(baseline_distro))
  df.dynamics <- data.frame(distance=c(dynamics_distro))

  # plot histograms of these distributions
  p.1 <- ggplot(df.baseline) +
    geom_histogram(aes(x=distance)) +
    xlab("baseline distance (Euclidean)")
  p.2 <- ggplot(df.dynamics) +
    geom_histogram(aes(x=distance)) +
    xlab("dynamics distance (Riemannian)")
  p <- grid.arrange(p.1, p.2, nrow=1)
  ggsave("output/distance_distributions.png", p, dpi=150, units="in", height=4, width=8)
}

n_hosts <- length(model_list$hosts)
unique_combos <- combn(1:n_hosts, 2)
n_unique_combos <- ncow(unique_combos)

quantile_mat.baseline <- matrix(NA, n_hosts, n_hosts)
quantile_mat.dynamics <- matrix(NA, n_hosts, n_hosts)

for(i in 1:n_hosts) {
  for(j in 1:n_hosts) {
    if(i < j) {
      # only consider unique combos
      distance <- baseline_dist[i,j]
      distance.quantile <- sum(baseline_distro < distance)/unique_combos
      quantile_mat.baseline[i,j] <- distance.quantile
      distance <- dynamics_dist[i,j]
      distance.quantile <- sum(dynamics_distro < distance)/unique_combos
      quantile_mat.dynamics[i,j] <- distance.quantile
    }
  }
}

# for random pairs
n_sample <- 1000
iterate_combos <- sample(1:n_unique_combos)[1:n_sample]
selected_dist.baseline <- c()
selected_dist.dynamics <- c()
for(i in 1:n_sample) {
  hosts <- c(unique_combos[1,iterate_combos[i]], unique_combos[2,iterate_combos[i]])
  selected_dist.baseline <- c(selected_dist.baseline, quantile_mat.baseline[min(hosts), max(hosts)])
  selected_dist.dynamics <- c(selected_dist.dynamics, quantile_mat.dynamics[min(hosts), max(hosts)])
}

# plot any evidence of correlation?
df <- data.frame(x=selected_dist.baseline, y=selected_dist.dynamics)

p <- ggplot(df) +
  geom_point(aes(x=x, y=y)) +
  xlab("baseline quantile") +
  ylab("dynamics quantile")
ggsave("output/distance_quantile_correlations.png", p, dpi=150, units="in", height=6, width=6)

# visualize via heatmaps
df <- gather_array(quantile_mat.baseline, "quantile", "host.1", "host.2")
p.1 <- ggplot(df, aes(x=host.1, y=host.2)) +
  geom_tile(aes(fill = quantile)) +
  scale_fill_viridis_c(option = "B", direction = -1) +
  xlab("host 1") +
  ylab("host 2") +
  theme_light()
df <- gather_array(quantile_mat.dynamics, "quantile", "host.1", "host.2")
p.2 <- ggplot(df, aes(x=host.1, y=host.2)) +
  geom_tile(aes(fill = quantile)) +
  scale_fill_viridis_c(option = "B", direction = -1) +
  xlab("host 1") +
  ylab("host 2") +
  theme_light()
p <- grid.arrange(p.1, p.2, nrow=2)
ggsave("output/quantile_heatmaps.png", p, dpi=150, units="in", height=8, width=5)

