library(ROL)
library(stringr)
library(driver)
library(phyloseq)
library(dplyr)

plot_shuffled <- function(a, b, lim1, lim2, max_pairs = NULL) {
  # Smaller indices roughly correspond to higher abundance taxa and larger indices to lower abundance taxa.
  # Also in the associations matrix, we remove redundant associations by reducing (1x2, 2x1) to (1x2) -- i.e.
  #   the first index will always be enriched for more abundant taxa and the second index for lowly abundant
  #   taxa. When plotting taxon pairs against each other, shuffle the first and second indices to remove this.
  # In other words, we'll randomly use a as x (b as y) for some indices and vice versa for others.
  if(!is.null(max_pairs)) {
    # Downsample
    selected <- sample(1:length(a))[1:max_pairs]
    a <- a[selected]
    b <- b[selected]
  }
  selected <- as.logical(sample(c(0, 1), size = length(a), replace = TRUE))
  x <- a
  y <- b
  x[selected] <- b[selected]
  y[selected] <- a[selected]
  plot(x, y, xlim = c(lim1, lim2), ylim = c(lim1, lim2))
}

tax_level <- "ASV"
Sigmas <- load_MAP_estimates(tax_level = tax_level, DLM = TRUE, logratio = "clr")
mat <- plot_interaction_heatmap(tax_level = tax_level, logratio = "clr", Sigmas = Sigmas, DLM = TRUE, cluster = TRUE, taxon_idx = NULL, show_plot = FALSE, return_matrix = TRUE)

str(mat, max.level = 1)

mean_assoc <- colMeans(mat$interactions.reordered)
plot(mean_assoc[order(mean_assoc)])

top_pos_idx <- which(mean_assoc > 0.2)
top_pos <- sapply(mat$labels.reordered[top_pos_idx], function(x) {
  as.numeric(str_match(x, "(\\d+)_(\\d+)")[,2:3])
})
colnames(top_pos) <- NULL

top_neg_idx <- which(mean_assoc < -0.2)
top_neg <- sapply(mat$labels.reordered[top_neg_idx], function(x) {
  as.numeric(str_match(x, "(\\d+)_(\\d+)")[,2:3])
})
colnames(top_neg) <- NULL

no_assoc <- which(abs(mean_assoc) < 0.05)
top_null <- sapply(mat$labels.reordered[no_assoc], function(x) {
  as.numeric(str_match(x, "(\\d+)_(\\d+)")[,2:3])
})
colnames(top_null) <- NULL

dim(top_pos)
dim(top_neg)
dim(top_null)

# (1) Is there evidence of one taxon wildly more prevalent than the others in terms of associations? (NO)
freqs_pos <- table(c(top_pos))
plot(as.numeric(names(freqs_pos)), as.numeric(freqs_pos), xlab = "taxon ID", ylab = "frequency", main = "positive correlator frequencies") # positive

freqs_neg <- table(c(top_neg))
plot(as.numeric(names(freqs_neg)), as.numeric(freqs_neg), xlab = "taxon ID", ylab = "frequency", main = "negative correlator frequencies") # negative

freqs_null <- table(c(top_null))
plot(as.numeric(names(freqs_null)), as.numeric(freqs_null), xlab = "taxon ID", ylab = "frequency", main = "null correlator frequencies") # null

# (2) Is there an obvious abundance pattern here?
host <- "ACA"
sample_fit <- readRDS(paste0("C:/Users/kim/Documents/ROL/output/model_fits/ASV_MAP/",host,"_labraduckfit.rds"))
counts <- sample_fit$Y # taxa x samples
clr.counts <- clr_array(counts + 0.5, parts = 1)
clr.avg <- rowMeans(clr.counts)

# For each pair of correlators, take plot the densities of the 
max_pairs <- ncol(top_pos)
shuffle_idx <- sample(1:ncol(top_null))[1:max_pairs]
df <- data.frame(log_abundance1 = clr.avg[top_null[1,shuffle_idx]], log_abundance2 = clr.avg[top_null[2,shuffle_idx]], type = 1)
shuffle_idx <- sample(1:ncol(top_pos))[1:max_pairs]
df <- rbind(df, data.frame(log_abundance1 = clr.avg[top_pos[1,shuffle_idx]], log_abundance2 = clr.avg[top_pos[2,shuffle_idx]], type = 2))
shuffle_idx <- sample(1:ncol(top_neg))[1:max_pairs]
df <- rbind(df, data.frame(log_abundance1 = clr.avg[top_neg[1,shuffle_idx]], log_abundance2 = clr.avg[top_neg[2,shuffle_idx]], type = 3))
df$type <- as.factor(df$type)
levels(df$type) <- c("null", "positive", "negative")
# shuffle column 1 & 2 within rows for the same reason noted in plot_shuffled()
for(i in 1:nrow(df)) {
  if(runif(1, min = 0, max = 1) > 0.5) {
    df[i,1:2] <- df[i,2:1]
  }
}
ggplot(df, aes(x = log_abundance1, y = log_abundance2)) +
  geom_point() +
  facet_wrap(~ type, ncol = 3)


# Positive correlators: Tend to lack things at different quantiles of expression (opposite sides of the mean).
# Negative correlators: Tend to be enriched for things at different quantiles of expression (same side).

# Pull taxonomy.
data <- load_data(tax_level = "ASV", host_sample_min = 75, count_threshold = 5, sample_threshold = 0.2)
params <- formalize_parameters(data)
tax <- get_taxonomy(data, alr_ref = params$alr_ref)

short_tax <- function(x) {
  name_idx <- max(which(!is.na(x)))
  paste0("ASV in ",names(x)[name_idx]," ",x[name_idx])
}

eval_seasonality <- function(eval_obj, data, host) {
  host <<- host
  # raw data
  host_data <- subset_samples(data, sname == host)
  host_md <- sample_data(host_data)
  # counts <- otu_table(host_data)@.Data
  # clr.counts <- clr(counts + 0.5) # samples x taxa
  host_data <- readRDS(paste0("output/model_fits/ASV_MAP/",host,"_labraduckfit.rds"))
  clr.fit <- to_clr(host_data$fit)
  clr.eta <- clr.fit$Eta[,,1]
  scores <- matrix(NA, ncol(eval_obj), 3)
  for(idx in 1:ncol(eval_obj)) {
    pair <- eval_obj[,idx]
    # scores[idx,] <- score(clr.counts[,pair[1]], clr.counts[,pair[2]], host_md$season)
    scores[idx,] <- score(clr.eta[pair[1],], clr.eta[pair[2],], host_md$season)
  }
  colMeans(scores)
}

# Quick and dirty check for "seasonality" using the score function from `evaluate_DLM_correlators.R`
data <- readRDS("input/filtered_ASV_5_20.rds")
# Filter to individuals with at least 75 samples.
md <- sample_data(data)
hosts <- suppressWarnings(unname(unlist(md %>%
                                          group_by(sname) %>%
                                          tally() %>%
                                          filter(n >= 75) %>%
                                          select(sname))))
data <- subset_samples(data, hosts %in% hosts)
md <- sample_data(data)

avg_scores_pos <- matrix(NA, length(hosts), 3)
for(h in 1:length(hosts)) {
  cat("Evaluating",hosts[h],"\n")
  avg_scores_pos[h,] <- eval_seasonality(top_pos, data, hosts[h])
}

avg_scores_neg <- matrix(NA, length(hosts), 3)
for(h in 1:length(hosts)) {
  cat("Evaluating",hosts[h],"\n")
  avg_scores_neg[h,] <- eval_seasonality(top_neg, data, hosts[h])
}

avg_scores_null <- matrix(NA, length(hosts), 3)
for(h in 1:length(hosts)) {
  cat("Evaluating",hosts[h],"\n")
  avg_scores_null[h,] <- eval_seasonality(top_null, data, hosts[h])
}

avg_scores <- avg_scores_null
rownames(avg_scores) <- hosts
df <- gather_array(avg_scores, "score", "host", "type")
df$type <- as.factor(df$type)
levels(df$type) <- c("agreement", "seasonality (taxon 1)", "seasonality (taxon 2)")

ggplot(df, aes(x = type, y = score)) +
  geom_boxplot() +
  ylim(c(0.3, 0.7))

df <- data.frame(score = c(avg_scores_null[,2], avg_scores_null[,3]), type = 1)
df <- rbind(df, data.frame(score = c(avg_scores_pos[,2], avg_scores_pos[,3]), type = 2))
df <- rbind(df, data.frame(score = c(avg_scores_neg[,2], avg_scores_neg[,3]), type = 3))
df$type <- as.factor(df$type)
levels(df$type) <- c("null", "positive", "negative")
ggplot(df, aes(x = type, y = score)) +
  geom_boxplot()

# positive_scores <- colMeans(avg_scores)
# negative_scores <- colMeans(avg_scores)
# null_scores <- colMeans(avg_scores)

positive_scores
negative_scores
null_scores
