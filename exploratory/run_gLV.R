library(ggplot2)
library(dplyr)
library(ROL)
library(matrixsampling)
library(gridExtra)

# Example scenarios
#' (0, 0.0): Nothing in common between hosts
#' (0, 0.5): 50% Environment in common
#' (1, 0.0): "Innate" parameters (baseline abundances and cooperative/competitive dynamics) in common but no environment
#' (2, 0.0): "Innate" parameters plus response to environment in common but no environment itself
#' (2, 0.5): "Innate" parameters plus response to environment + 50% environment in common

S <- 10
H <- 10
reps <- 20

save_file <- "simulation_sweep_5.png"
almean <- 0.8
shared_param_scenarios <- c(0, 1, 2)
shared_env_scenarios <- c(0.0, 0.5, 0.6, 0.7, 0.8, 0.9)

p1 <- almean
df <- data.frame(x = c(), y = c(), shared_param_scenario = c(), percent_shared_environment = c(), replicate = c())
for(p2 in shared_param_scenarios) {
  for(p3 in shared_env_scenarios) {
    cat("Generating replicates for scenario: (",p1,",",p2,",",p3,")\n")
    for(j in 1:reps) {
      cat("\tReplicate:",j,"\n")
      sim <- simulate_scenario(S = S, H = H,
                               almean = p1,
                               shared_param_level = p2,
                               shared_noise_proportion = p3)
      summary_xy <- summarize_all_pairs_2D(sim$heatmap)
      df <- rbind(df,
                  data.frame(x = summary_xy$x,
                             y = summary_xy$y,
                             shared_param_scenario = p2,
                             percent_shared_environment = p3,
                             replicate = j))
    }
  }
}

df$shared_param_scenario <- as.factor(df$shared_param_scenario)
df$percent_shared_environment <- as.factor(df$percent_shared_environment)
  
palette <- generate_highcontrast_palette(length(levels(df$percent_shared_environment)))
plot_list <- list()
for(i in 0:2) {
  p <- ggplot(df[df$shared_param_scenario == i,], aes(x = x, y = y, color = percent_shared_environment)) +
    geom_point(size = 2) +
    xlim(-0.1, 1.1) +
    ylim(0, 1) +
    theme(legend.position = "bottom") +
    scale_color_manual(values = palette)
  plot_list[[i+1]] <- p
}
p <- grid.arrange(grobs = plot_list, ncol = 3)
ggsave(save_file, p, units = "in", dpi = 100, height = 4, width = 11)

# plot_gLV_bars(sim$series[[1]])
# plot_gLV_heatmap(sim$heatmap)
# plot_gLV_hockeystick(sim$heatmap)









