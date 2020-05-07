library(ROL)
library(ggplot2)
library(dplyr)
library(purrr)
library(driver)
library(stringr)
library(phyloseq)
library(mixtools)

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
if(FALSE) {
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

  # gotchas
  x <- c(0.05, 0.04, 0.03, 0.8, 0.9, 1.0)
  cat("Score:",round(calc_universality_score(x), 3),"\n")
  
  x <- c(-0.05, -0.04, -0.03, 0.8, 0.9, 1.0)
  cat("Score:",round(calc_universality_score(x), 3),"\n")
}

# visualize scores against "The Rug"
if(FALSE) {
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
}

# visualize a null distribution of scores
if(FALSE) {
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
  ggsave(file.path(save_dir,"null_distribution_universality.png"), p, units = "in", dpi = 100, height = 5, width = 8) 
}

# calculate standard deviation part of universality score
calc_universality_score.1 <- function(x) {
  max_sd <- sd(c(rep(-1, length(x)/2), rep(1, length(x)/2)))
  return(((max_sd - sd(x))/max_sd))
}

# calculate effect size part of universality score
calc_universality_score.2 <- function(x) {
  return(abs(mean(x)))
}

# aggregate these to per-microbe scores and plot distributions
hosts <- names(Sigmas)
Sigmas_corr <- Sigmas
for(host in hosts) {
  Sigmas_corr[[host]] <- cov2cor(Sigmas_corr[[host]])
}

all_scores <- c()
# all_scores.1 <- c()
# all_scores.2 <- c()
score_labels <- c()
median_scores <- c()
for(focal_microbe in 1:n_taxa) {
  focal_scores <- c()
  # focal_scores.1 <- c()
  # focal_scores.2 <- c()
  for(other_microbe in setdiff(1:n_taxa, focal_microbe)) {
    focal_interactions <- c()
    for(host in hosts) {
      focal_interactions <- c(focal_interactions, Sigmas_corr[[host]][focal_microbe, other_microbe])
    }
    focal_scores <- c(focal_scores, calc_universality_score(focal_interactions))
    # focal_scores.1 <- c(focal_scores.1, calc_universality_score.1(focal_interactions))
    # focal_scores.2 <- c(focal_scores.2, calc_universality_score.2(focal_interactions))
  }
  all_scores <- c(all_scores, focal_scores)
  # all_scores.1 <- c(all_scores.1, focal_scores.1)
  # all_scores.2 <- c(all_scores.2, focal_scores.2)
  score_labels <- c(score_labels, rep(focal_microbe, length(focal_scores)))
  median_scores[focal_microbe] <- median(focal_scores)
}

# rank all by universality
data <- load_data(tax_level=tax_level)
log_abundances <- log(otu_table(data)@.Data + 0.5)
avg_log_abundance <- colMeans(log_abundances)
avg_log_deviance <- apply(log_abundances, 2, sd)
names(avg_log_abundance) <- NULL
# reorder the taxa to match the ordering in the model output
alr_ref <- formalize_parameters(data)$alr_ref
avg_log_abundance <- avg_log_abundance[c(setdiff(1:n_taxa, alr_ref), alr_ref)]
avg_log_deviance <- avg_log_deviance[c(setdiff(1:n_taxa, alr_ref), alr_ref)]
taxonomy <- get_taxonomy(data, alr_ref)

reorder <- order(median_scores, decreasing = TRUE)
sink(file.path(save_dir,"microbe_universality_table.txt"))
cat("Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tScore\tAvg. log abundance\tStd. dev. log abundance\n")
for(rr in reorder) {
  cat(paste0(taxonomy[rr,1],"\t",taxonomy[rr,2],"\t",taxonomy[rr,3],"\t",taxonomy[rr,4],"\t",taxonomy[rr,5],"\t",taxonomy[rr,6],"\t"))
  cat(paste0(round(median_scores[rr], 3),"\t"))
  cat(paste0(round(avg_log_abundance[rr], 3),"\t"))
  cat(paste0(round(avg_log_deviance[rr], 3),"\n"))
}
sink()

# full score
df <- data.frame(score = all_scores, focal_microbe = score_labels)
p <- ggplot(df) +
  geom_density(aes(x = score)) +
  xlim(0, 0.7) +
  facet_wrap(vars(focal_microbe), nrow = 6)
ggsave(file.path(save_dir, "microbe_universality_scores.png"), p, units = "in", dpi = 100, height = 10, width = 20)

# # part 1
# df.1 <- data.frame(score = all_scores.1, focal_microbe = score_labels)
# p <- ggplot(df.1) +
#   geom_density(aes(x = score)) +
#   #xlim(0, 0.7) +
#   facet_wrap(vars(focal_microbe), nrow = 6)
# ggsave(file.path(save_dir, "microbe_universality_scores_pt1.png"), p, units = "in", dpi = 100, height = 10, width = 20)

# # part 2
# df.2 <- data.frame(score = all_scores.2, focal_microbe = score_labels)
# p <- ggplot(df.2) +
#   geom_density(aes(x = score)) +
#   #xlim(0, 0.7) +
#   facet_wrap(vars(focal_microbe), nrow = 6)
# ggsave(file.path(save_dir, "microbe_universality_scores_pt2.png"), p, units = "in", dpi = 100, height = 10, width = 20)

# visualize the rug with these scores
tax <- assign_concise_taxonomy(tax_level = tax_level, logratio = logratio)
most_universal_idx <- which(median_scores == max(median_scores))
cat("Max universal taxon:",tax[most_universal_idx],", index:",most_universal_idx,"\n")
plot_interaction_heatmap(tax_level = tax_level, logratio = logratio, Sigmas = Sigmas, taxon_idx = most_universal_idx)
med_universal_idx <- which(order(median_scores) == 50)
cat("Med universal taxon:",tax[med_universal_idx],", index:",med_universal_idx,"\n")
plot_interaction_heatmap(tax_level = tax_level, logratio = logratio, Sigmas = Sigmas, taxon_idx = med_universal_idx)
least_universal_idx <- which(median_scores == min(median_scores))
cat("Min universal taxon:",tax[least_universal_idx],", index:",least_universal_idx,"\n")
plot_interaction_heatmap(tax_level = tax_level, logratio = logratio, Sigmas = Sigmas, taxon_idx = least_universal_idx)

# BIMODALITY
pairs_obj <- get_pairwise_correlations(tax_level=tax_level, logratio="clr")
labels <- pairs_obj$labels
interactions <- pairs_obj$interactions

# fit a 1 and 2-component 1D Gaussian mixture and evaluate evidence for bimodality viw a LRT test
# this takes ~1 min to run
significance <- rep(NA, ncol(interactions))
for(i in 1:ncol(interactions)) {
  cat(i,"\n")
  x <- interactions[,i]
  result = tryCatch({
    res.1 <- sum(dnorm(x, mean = mean(x), sd = sd(x), log = TRUE))
    res.2 <- sum(normalmixEM(x, k = 2)$loglik)
    lrt_val <- -2*(res.1 - res.2)
    significance[i] <- pchisq(lrt_val, df = 3, lower.tail=FALSE)
  },
  error=function(cond) {
    significance[i] <- NA
  },
  warning=function(cond) {
    significance[i] <- NA
  })
}

# order interactions by evidence of bimodality
bimodal_evidence <- order(significance)

# evidence of bimodality in interactions?
df <- gather_array(interactions, "correlation", "host", "interaction")
p <- ggplot(df[df$interaction %in% bimodal_evidence[1:200],]) +
  geom_density(aes(x = correlation)) +
  facet_wrap(vars(interaction), nrow = 10)
ggsave("test.png", p, units = "in", dpi = 100, height = 10, width = 20)

# these were interesting "bimodal" interactions identified by hand
interesting_pair <- labels[5551] # also 3073, 3358
microbe_pair <- as.numeric(strsplit(interesting_pair, "_")[[1]])

# get some individuals centered around the two modes
correlations <- interactions[,5551]
modes <- normalmixEM(correlations, k = 2)$mu
host_mode1 <- which.min(abs(correlations - modes[1]))
host_mode2 <- which.min(abs(correlations - modes[2]))

hosts <- names(Sigmas)
plot_posterior_predictive(host=hosts[host_mode1], tax_level="ASV", predict_coords=c(microbe_pair[1], microbe_pair[2]))
plot_posterior_predictive(host=hosts[host_mode2], tax_level="ASV", predict_coords=c(microbe_pair[1], microbe_pair[2]))
