library(ROL)
library(ggplot2)
library(dplyr)
library(purrr)
library(driver)

tax_level <- "ASV"
logratio <- "clr"

uncertain_scores <- FALSE

if(uncertain_scores) {
  Sigmas <- load_full_posteriors(tax_level="ASV", logratio="clr")
  n_hosts <- length(Sigmas)
  n_taxa <- dim(Sigmas[[1]])[1]
  n_interactions <- (n_taxa^2)/2 - n_taxa/2
  n_samples <- dim(Sigmas[[1]])[3]
  interactions <- array(NA, dim=c(n_hosts, n_interactions, n_samples))
  host_reordering <- NULL
  interaction_reordering <- NULL
  for(k in 1:n_samples) {
    cat("Sample",k,"\n")
    partial_Sigmas <- Sigmas
    for(i in 1:n_hosts) {
      partial_Sigmas[[i]] <- partial_Sigmas[[i]][,,k]
    }
    interactions[,,k] <- plot_interaction_heatmap(tax_level=tax_level, logratio=logratio, Sigmas=partial_Sigmas,
                                                  cluster=FALSE, return_matrix=TRUE)
    # calculate an ordering based on the first sample
    if(k == 1) {
      cat("Calculating an ordering...\n")
      d <- dist(interactions[,,k])
      clustering.hosts <- hclust(d)
      d <- dist(t(interactions[,,k]))
      clustering.interactions <- hclust(d)
      # reorder all
      host_reordering <- clustering.hosts$order
      interaction_reordering <- clustering.interactions$order
    }
    interactions[,,k] <- interactions[host_reordering,interaction_reordering,k]
  }
} else {
  Sigmas <- load_MAP_estimates(tax_level="ASV", logratio="clr")
  n_hosts <- length(Sigmas)
  n_taxa <- dim(Sigmas[[1]])[1]
  n_interactions <- (n_taxa^2)/2 - n_taxa/2
  interactions <- plot_interaction_heatmap(tax_level=tax_level, logratio=logratio, Sigmas=Sigmas,
                                           cluster=FALSE, return_matrix=TRUE)
  cat("Calculating an ordering...\n")
  d <- dist(interactions)
  clustering.hosts <- hclust(d)
  d <- dist(t(interactions))
  clustering.interactions <- hclust(d)
  # reorder all
  host_reordering <- clustering.hosts$order
  interaction_reordering <- clustering.interactions$order
  interactions <- interactions[host_reordering,interaction_reordering]
}

if(uncertain_scores) {
  save_dir <- check_output_dir(c("output","plots",tax_level))
} else {
  save_dir <- check_output_dir(c("output","plots",paste0(tax_level,"_MAP")))
}

# sanity checking with some exemplar scores
if(TRUE) {
  # perfect agreement
  x <- rep(1, 6)
  cat("Score:",round(calc_universality_score(x), 2),"\n")

  # perfect disagreement
  x <- c(1, 1, 1, -1, -1, -1)
  cat("Score:",round(calc_universality_score(x), 2),"\n")

  # high mean, low-ish variance
  x <- c(1, 0.95, 0.9, 0.85, 0.8, 0.75)
  cat("Score:",round(calc_universality_score(x), 2),"\n")

  # near-zero mean, low-ish variance
  x <- c(-0.15, -0.1, -0.05, 0, 0.05, 0.01)
  cat("Score:",round(calc_universality_score(x), 2),"\n")

  # near-zero mean, low-ish variance
  x <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)
  cat("Score:",round(calc_universality_score(x), 2),"\n")
}

quit()

if(uncertain_scores) {
  # plot heatmap/"rug" for the first few posterior samples
  for(k in 1:4) {
    df <- gather_array(interactions[,,1], "correlation", "host", "interaction")
    p <- ggplot(df, aes(interaction, host)) +
      geom_tile(aes(fill = correlation)) +
      scale_fill_gradient2(low = "darkblue", high = "darkred")
    ggsave(file.path(save_dir,paste0("universal_interactions_posterior_sample_",k,".png")),
           p, units="in", dpi=100, height=5, width=15)
  }

  # calculate intervals of uncertainty associated with scores
  # the intervals will come from posterior samples
  uncertain_scores_df <- data.frame(score=c(), interaction_idx=c(), sample=c())
  for(k in 1:n_samples) {
      uncertain_scores_df <- rbind(uncertain_scores_df, 
                                    data.frame(score=apply(interactions[,,k], 2, calc_universality_score),
                                              interaction_idx=1:n_interactions, sample=k))
  }
  
  # get 95% interval
  post_quantiles <- uncertain_scores_df %>%
    group_by(interaction_idx) %>%
    summarise(p2.5 = quantile(score, prob=0.025),
              mean = mean(score),
              p97.5 = quantile(score, prob=0.975)) %>%
    ungroup()
  post_quantiles <- as.data.frame(post_quantiles)
  
  # plot scores
  p <- ggplot(post_quantiles) +
    geom_ribbon(aes(x=interaction_idx, ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
    geom_line(aes(x=interaction_idx, y=mean), color="blue", group=1, size=0.25) +
    theme_minimal()
  ggsave(file.path(save_dir,paste0("universal_interactions_posterior_scores.png")),
          p, units="in", dpi=100, height=2, width=15)
} else {
  # plot MAP heatmap/"rug"
  df <- gather_array(interactions, "correlation", "host", "interaction")
  p <- ggplot(df, aes(interaction, host)) +
    geom_tile(aes(fill = correlation)) +
    scale_fill_gradient2(low = "darkblue", high = "darkred")
  ggsave(file.path(save_dir,paste0("universal_interactions_MAP.png")),
          p, units="in", dpi=100, height=5, width=15)

  # plot scores
  scores <- apply(interactions, 2, calc_universality_score)
  df <- data.frame(interaction=1:n_interactions, score=scores)
  p <- ggplot(df) +
    geom_line(aes(x=interaction, y=score), color="blue", group=1, size=0.25) +
    theme_minimal()
  ggsave(file.path(save_dir,paste0("universal_interactions_MAP_scores.png")),
          p, units="in", dpi=100, height=2, width=15)
}

if(FALSE) {
  # null distribution of scores?
  # permute interactions across hosts within a column
  permuted_interactions <- interactions
  n_interactions <- ncol(interactions)
  for(i in 1:n_hosts) {
    shuffle_order <- sample(1:n_interactions)
    permuted_interactions[i,] <- permuted_interactions[i,shuffle_order]
  }
  H0_scores <- c()
  for(j in 1:n_interactions) {
    H0_scores <- c(H0_scores, calc_universality_score(permuted_interactions[,j]))
  }

  df <- data.frame(x = H0_scores, which = "null")
  df <- rbind(df, data.frame(x = apply(interactions, 2, calc_universality_score), which = "actual"))
  p <- ggplot(df) +
    geom_density(aes(x = x, color = which)) +
    xlim(c(0, 1))
  ggsave(file.path(save_dir,"hull_distribution_universality.png"), p, units = "in", dpi = 100, height = 5, width = 8) 
}
