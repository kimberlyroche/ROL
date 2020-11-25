library(ROL)

# Scenarios (copied from generate_sim_params)
#' 1: Nothing in common between hosts
#' 2: 50% Environment in common
#' 3: 80% Environment in common
#' 4: "Innate" parameters (baseline abundances and cooperative/competitive dynamics) in common
#' 5: "Innate" parameters plus response to environment in common
#' 6: "Innate" parameters plus response to environment + 50% environment in common
#' 7: "Innate" parameters plus response to environment + 80% environment in commontemp <- generate_sim_params(7)

S <- 10
H <- 10
scenario <- 1
sim <- simulate_scenario(S = S, H = H, scenario = scenario)

plot_gLV_heatmap(sim$heatmap)
plot_gLV_hockeystick(sim$heatmap)
plot_gLV_bars(sim$series[[1]])
