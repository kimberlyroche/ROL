## ------------------------------------------------------------------------------------------------
##   VISUALIZING THE UNIVERSALITY SCORE
## ------------------------------------------------------------------------------------------------

library(ggplot2)
library(stringr)
library(ROL)
library(dplyr)
library(LaplacesDemon)
library(driver)
library(matrixsampling)
library(emdist)

calc_xy <- function(vSigmas) {
  if(is.null(dim(vSigmas))) {
    dim(vSigmas) <- c(length(vSigmas), 1)
  }
  shared_direction <- mean(apply(sign(vSigmas), 2, function(z) {
    max(table(z)/length(z))
  }))
  mean_strength <- mean(abs(apply(vSigmas, 2, mean)))
  list(x = shared_direction, y = mean_strength)
}

generate_counts <- function(num_taxa, num_samples, num_hosts, scenario, population_baseline_composition = NULL) {
  synthetic_counts <- array(NA, dim = c(num_taxa, num_samples, num_hosts))
  for(h in 1:num_hosts) {
    if(scenario == 1) {
      # total randomness; no individual baseline
      props_h <- rdirichlet(n = num_samples, alpha = rep(0.1, num_taxa))
      counts_h <- t(apply(props_h, 1, function(x) {
        rmultinom(1, size = rpois(1, 10000), prob = x)
      })) # dimensions are samples x taxa
    }
    if(scenario == 2) {
      # a host baseline gives correlation between samples
      baseline_h <- rdirichlet(n = 1, alpha = rep(1, num_taxa))
      counts_h <- matrix(NA, num_taxa, num_samples)
      for(n in 1:num_samples) {
        temp <- rdirichlet(n = 1, alpha = baseline_h*200)
        counts_h[,n] <- rmultinom(1, size = rpois(1, 10000), prob = temp)[,1]
      }
    }
    if(scenario == 3) {
      # population and host baselines correlate samples now
      # baseline_h <- rdirichlet(n = 1, alpha = population_baseline_composition*200)
      counts_h <- matrix(NA, num_taxa, num_samples)
      for(n in 1:num_samples) {
        # temp <- rdirichlet(n = 1, alpha = baseline_h*50)
        temp <- rdirichlet(n = 1, alpha = population_baseline_composition*50)
        counts_h[,n] <- rmultinom(1, size = rpois(1, 10000), prob = temp)[,1]
      }
    }
    synthetic_counts[,,h] <- counts_h
  }
  synthetic_counts
}

# host_contrib = 0.5 gives a 50% population-level, 50% host-level covariance structure to logratio taxa
generate_counts_w_cov <- function(num_taxa, num_samples, num_hosts, host_contrib = 0.5) {
  synthetic_counts <- array(NA, dim = c(num_taxa, num_samples, num_hosts))
  simulated_covariance_pop <- rinvwishart(1, num_taxa + 10, diag(num_taxa)*10)[,,1]
  for(h in 1:num_hosts) {
    # host level correlation
    simulated_covariance_host <- rinvwishart(1, num_taxa + 10, diag(num_taxa)*10)[,,1]
    simulated_covariance <- (1 - host_contrib)*simulated_covariance_pop + host_contrib*simulated_covariance_host
    counts_h <- matrix(NA, num_taxa, num_samples)
    logratios <- rmatrixnormal(1, matrix(0, num_taxa, num_samples), simulated_covariance, diag(num_samples))[,,1]
    for(n in 1:num_samples) {
      counts_h[,n] <- rmultinom(1, size = rpois(1, 10000), prob = clrInv(logratios[,n]))[,1]
    }
    synthetic_counts[,,h] <- counts_h
  }
  synthetic_counts
}

## ------------------------------------------------------------------------------------------------
##   Plot (empty) heatmap with scores
## ------------------------------------------------------------------------------------------------

shared_direction_proportion <- seq(from = 0.5, to = 1, by = 0.05)
mean_association_strength <- seq(from = 0, to = 1, by = 0.05)

data <- data.frame(x = c(), y = c(), z = c())

for(x in shared_direction_proportion) {
  for(y in mean_association_strength) {
    data <- rbind(data, data.frame(x = x, y = y, z = x*y))
  }
}

ggplot() +
  geom_raster(data = data, aes(x = x, y = y, fill = z)) +
  geom_text(data = data, aes(x = x, y = y, label = format(round(z, 2), nsmall = 2, color = "#333333"))) +
  scale_fill_gradient(low = "white", high = "red") +
  # geom_rect(aes(xmin = 0.775, xmax = 0.825, ymin = 0.125, ymax = 0.175), color = "black", fill = NA, alpha = 1) +
  xlab("proportion shared direction") +
  ylab("mean strength of associations") +
  labs(fill = "weighted \nuniversality\nscore")

## ------------------------------------------------------------------------------------------------
##   Plot ABRP data set
## ------------------------------------------------------------------------------------------------

Sigmas <- load_MAP_estimates(tax_level = "ASV", logratio = "clr")
Sigmas <- lapply(Sigmas, function(x) cov2cor(x))

n_hosts <- length(Sigmas)
P <- dim(Sigmas[[1]])[1]

vectorized_Sigmas <- matrix(NA, n_hosts, (P^2)/2 - P/2)
for(i in 1:n_hosts){
  vectorized_Sigmas[i,] <- Sigmas[[i]][upper.tri(Sigmas[[i]], diag = F)]
}

interactions <- get_pairwise_correlations(tax_level = "ASV", logratio = "clr", Sigmas = Sigmas)

selected_interactions <- list()
for(i in 1:nrow(Sigmas[[1]])) {
  selected_interactions[[i]] <- which(str_detect(interactions$labels, regex(paste0("(^",i,"_.*$|^.*?_",i,"$)"), ignore_case = TRUE)))
}

point_data <- data.frame(x = c(), y = c(), taxon = c())
for(taxon in 1:length(selected_interactions)) {
  cat("Taxon:",taxon,"\n")
  for(i in selected_interactions[[taxon]]) {
    temp <- calc_xy(vectorized_Sigmas[,i])
    point_data <- rbind(point_data, data.frame(x = temp$x, y = temp$y, taxon = taxon))
  }
}

temp <- point_data %>%
  group_by(taxon) %>%
  summarise(mean = mean(y)) %>%
  arrange(desc(mean))

point_data$selected <- FALSE
point_data[point_data$taxon == temp[1,]$taxon,]$selected <- TRUE # highest avg. mean association strength
point_data[point_data$taxon == temp[nrow(temp),]$taxon,]$selected <- TRUE # lowest ...

abrp_data <- load_data(tax_level = "ASV")
alr_ref <- formalize_parameters(abrp_data)$alr_ref
taxonomy <- get_taxonomy(abrp_data, alr_ref)
taxonomy[temp[1,]$taxon,]
taxonomy[temp[nrow(temp),]$taxon,]

ggplot() +
  geom_raster(data = data, aes(x = x, y = y, fill = z)) +
  scale_fill_gradient(low = "white", high = "red") +
  geom_point(data = point_data[point_data$selected == FALSE,], aes(x = x, y = y), color = "#444444") +
  geom_point(data = point_data[point_data$selected == TRUE,], aes(x = x, y = y), color = "#16a3db", size = 3) +
  xlab("proportion shared direction") +
  ylab("mean strength of associations") +
  labs(fill = "weighted \nuniversality\nscore")

## ------------------------------------------------------------------------------------------------
##   Simulations
## ------------------------------------------------------------------------------------------------

# Functions

# Plot the 2D heatmap of universality scores and map pairwise interactions onto these as a point
# cloud
visualize_synthetic_data <- function(synthetic_point_data, show_rug = FALSE) {
  # plot heatmap with points
  p <- ggplot() +
    geom_raster(data = data, aes(x = x, y = y, fill = z)) +
    scale_fill_gradient(low = "white", high = "red") +
    geom_point(data = synthetic_point_data, aes(x = x, y = y), color = "#444444") +
    xlab("proportion shared direction") +
    ylab("mean strength of associations") +
    labs(fill = "weighted \nuniversality\nscore")
  show(p)
  
  if(show_rug) {
    # render the "rug" for this case to get another view of what's going on; takes ~30 sec. to run
    d <- dist(synthetic_Sigmas)
    clustering.hosts <- hclust(d)
    d <- dist(t(synthetic_Sigmas))
    clustering.interactions <- hclust(d)
    interactions.reordered <- synthetic_Sigmas[clustering.hosts$order,]
    interactions.reordered <- interactions.reordered[,clustering.interactions$order]
    
    plot_df <- gather_array(interactions.reordered, "correlation", "host", "interaction")
    p <- ggplot(plot_df, aes(x = interaction, y = host, fill = correlation)) +
      geom_tile() +
      scale_fill_gradient2(low = "darkblue", high = "darkred")
    show(p)
  }
}

# Calculate the 2D "universality" mapping coordinates for all simulated interactions
map_synthetic_data <- function(counts) {
  D <- dim(counts)[1]
  N_h <- dim(counts)[2]
  H <- dim(counts)[3]
  synthetic_Sigmas <- matrix(NA, H, (D^2/2) - D/2)
  # now calculate correlation across pairs
  # at D = 150, H = 50, N_h = 100 this takes about 1 min.
  for(h in 1:H) {
    clr_h <- clr(t(counts[,,h]) + 0.5)
    corr_h <- c()
    for(i in 1:(D_new - 1)) {
      for(j in (i+1):D_new) {
        temp <- cor(clr_h[,i], clr_h[,j])
        corr_h <- c(corr_h, temp)
      }
    }
    synthetic_Sigmas[h,] <- corr_h
  }
  synthetic_point_data <- data.frame(x = c(), y = c())
  for(i in 1:ncol(synthetic_Sigmas)) {
    # if(i %% 100 == 0) {
    #   cat("Interaction:",i,"\n")
    # }
    temp <- calc_xy(synthetic_Sigmas[,i])
    synthetic_point_data <- rbind(synthetic_point_data, data.frame(x = temp$x, y = temp$y))
  }
  synthetic_point_data
}

# Calculate Earth Mover's Distance between synthetic and observed interaction point clouds
calculate_distance <- function(synthetic_df, observed_df) {
  # bin the 2D map of universality
  # we'll use Earth Mover's Distance to assess how similar two point cloud distributions are
  binned_x <- unname(table(cut(synthetic_df$x, breaks = seq(0.5, 1, length.out = 20))))
  binned_y <- unname(table(cut(synthetic_df$y, breaks = seq(0, 1, length.out = 20))))
  synthetic_binned_data <- rbind(binned_x, binned_y)
  
  binned_x <- unname(table(cut(observed_df$x, breaks = seq(0.5, 1, length.out = 20))))
  binned_y <- unname(table(cut(observed_df$y, breaks = seq(0, 1, length.out = 20))))
  observed_binned_data <- rbind(binned_x, binned_y)
  
  emd2d(synthetic_binned_data, observed_binned_data)
}

# visualize the real data
visualize_synthetic_data(point_data)
ggsave("real_data.png", units = "in", dpi = 150, height = 7, width = 9)

D <- dim(Sigmas[[1]])[1]
D <- 150
H <- length(Sigmas)
H <- 50
N_h <- 100

# scenario <- 1
# 1: all compositions are random within and between hosts
#    taxa are independent
# 2: hosts have a baseline composition
#    taxa are independent
# 3: the population has a baseline composition and hosts have baseline compositions correlated with this
#    taxa are independent
# baseline_p <- rdirichlet(n = 1, alpha = rep(1, D))
# counts <- generate_counts(D, N_h, H, scenario = scenario, population_baseline_composition = baseline_p)

# 4: some combination of population- and host-level correlated taxa
sweep_host_contrib <- seq(from = 0, to = 1, by = 0.1)
for(host_contrib in sweep_host_contrib) {
  counts <- generate_counts_w_cov(D, N_h, H, host_contrib = host_contrib)
  synthetic_point_data <- map_synthetic_data(counts)
  # visualize_synthetic_data(synthetic_point_data)
  dist_value <- calculate_distance(synthetic_point_data, point_data)
  cat("Distance for combo:",round(1 - host_contrib, 1),"x",round(host_contrib, 1),":",round(dist_value, 2),"\n")
  visualize_synthetic_data(synthetic_point_data)
  ggsave(paste0("simulated_data_",host_contrib,".png"), units = "in", dpi = 150, height = 7, width = 9)
}


