library(ggplot2)
library(dplyr)
library(ROL)

# Scenarios
#' 1: Nothing in common between hosts
#' 2: 50% Environment in common
#' 3: 80% Environment in common
#' 4: "Innate" parameters (baseline abundances and cooperative/competitive dynamics) in common
#' 5: "Innate" parameters plus response to environment in common
#' 6: "Innate" parameters plus response to environment + 50% environment in common
#' 7: "Innate" parameters plus response to environment + 60% environment in common
#' 8: "Innate" parameters plus response to environment + 70% environment in common
#' 9: "Innate" parameters plus response to environment + 80% environment in common

S <- 10
H <- 10
reps <- 20
param_scenarios <- data.frame(p1 = rep(2, 6), p2 = c(0, 0.5, 0.6, 0.7, 0.8, 0.9))

df <- data.frame(x = c(), y = c(), scenario = c(), replicate = c())

for(i in 1:nrow(param_scenarios)) {
  cat("Generating replicates for scenario",i,"\n")
  for(j in 1:reps) {
    sim <- simulate_scenario(S = S, H = H,
                             shared_param_level = param_scenarios[i,]$p1,
                             shared_noise_proportion = param_scenarios[i,]$p2)
    summary_xy <- summarize_all_pairs_2D(sim$heatmap)
    df <- rbind(df, data.frame(x = summary_xy$x, y = summary_xy$y, scenario = i, replicate = j))
  }
}
df$scenario <- as.factor(df$scenario)
palette <- generate_highcontrast_palette(nrow(param_scenarios))
p <- ggplot(df, aes(x = x, y = y, color = scenario)) +
  geom_point(size = 2) +
  xlim(-0.1, 1.1) +
  ylim(0, 1) +
  scale_color_manual(values = palette)
show(p)

# plot_gLV_bars(sim$series[[1]])
# plot_gLV_heatmap(sim$heatmap)
# plot_gLV_hockeystick(sim$heatmap)

