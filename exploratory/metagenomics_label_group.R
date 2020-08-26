library(phyloseq)
library(ggplot2)

selected_hosts <- c("DUN", "EAG", "HON", "LAZ", "NOB", "ONY", "ORI", "VIN", "VOT")
data <- readRDS("input/ps0.rds")

plot_df <- data.frame(collection_date = c(), host = c(), group = c())
for(host in selected_hosts) {
  subset_data <- subset_samples(data, sname == host)
  subset_md <- sample_data(subset_data)
  plot_df <- rbind(plot_df, data.frame(collection_date = subset_md$collection_date,
                                       host = host,
                                       social_group = subset_md$grp))
}

# plot_df$social_group <- sapply(plot_df$social_group, function(x) substr(x, 1, 3))
plot_df$host <- as.factor(plot_df$host)
plot_df$social_group <- as.factor(plot_df$social_group)
plot_df$collection_date <- as.Date(plot_df$collection_date)
head(plot_df)

cbPalette <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(plot_df, aes(x = collection_date, y = host, color = social_group)) +
  geom_point(size = 3) +
  scale_colour_manual(values=cbPalette)
