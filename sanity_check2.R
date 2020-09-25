# This is a sanity check per Lawrence's suggestion. Are the patterns of correlation we see in the
# "rug" driven by compositional changes?

library(ROL)
library(ggplot2)
library(driver)
library(phyloseq)
library(stringr)

tax_level <- "ASV"
logratio <- "clr"
model_list <- get_fitted_model_list(tax_level = tax_level, MAP = TRUE)
Sigmas <- load_MAP_estimates(tax_level = tax_level, logratio = "clr")
pairs_obj <- get_pairwise_correlations(tax_level = tax_level, logratio = logratio, Sigmas = Sigmas)
labels <- pairs_obj$labels
interactions <- pairs_obj$interactions

D <- dim(Sigmas[[1]])[1]
H <- length(Sigmas)
interactions.p <- matrix(NA, dim(interactions)[1], dim(interactions)[2])
for(m in 1:H) {
  counts <- readRDS(model_list$model_list[m])$Y
  proportions <- apply(counts, 2, function(x) x/sum(x))
  it <- 1
  for(i in 2:D) { # iterate columns
    for(j in 1:(i-1)) { # iterate rows
      # taxon i, taxon j
      interactions.p[m,it] <- cor(proportions[j,], proportions[i,])
      it <- it + 1
    }
  }
  # interactions.p is in the same order as interactions
}

d <- dist(interactions)
clustering.hosts <- hclust(d)
d <- dist(t(interactions))
clustering.interactions <- hclust(d)
# reorder all
interactions.reordered <- interactions[clustering.hosts$order,]
interactions.reordered <- interactions.reordered[,clustering.interactions$order]
labels.reordered <- labels[clustering.interactions$order]
interactions.p.reordered <- interactions.p[clustering.hosts$order,]
interactions.p.reordered <- interactions.p.reordered[,clustering.interactions$order]

df <- gather_array(interactions.reordered, "correlation", "host", "pair")
p <- ggplot(df, aes(pair, host)) +
  geom_tile(aes(fill = correlation)) +
  scale_fill_gradient2(low = "darkblue", high = "darkred")
ggsave(paste0("microbe_pair_correlations_1.png"), p, units="in", dpi=150, height=5, width=15)

df <- gather_array(interactions.p.reordered, "correlation", "host", "pair")
p <- ggplot(df, aes(pair, host)) +
  geom_tile(aes(fill = correlation)) +
  scale_fill_gradient2(low = "darkblue", high = "darkred")
ggsave(paste0("microbe_pair_correlations_2.png"), p, units="in", dpi=150, height=5, width=15)

# In modeled data, find top 100 (+) correlators and top 100 (-) correlators.
pairs <- ncol(interactions.reordered)
mean_associations <- colMeans(interactions.reordered)
mean_associations.ordered <- order(mean_associations)
top_neg <- mean_associations.ordered[1:50]
top_pos <- mean_associations.ordered[(pairs-49):pairs]
# mean_associations[top_neg]
# mean_associations[top_pos]

data <- load_data(tax_level = tax_level)
counts <- t(otu_table(data)@.Data)
mean_log_abundances <- unname(apply(counts, 1, function(x) mean(log(x + 0.5))))
quartiles <- cut(mean_log_abundances, c(0, quantile(mean_log_abundances, probs = c(0.25, 0.5, 0.75, 1))))

lower_pair <- c()
upper_pair <- c()
for(tax_pair in top_neg) {
  num_pair <- as.numeric(str_split(labels.reordered[tax_pair], "_")[[1]])
  q1 <- as.numeric(quartiles[num_pair[1]])
  q2 <- as.numeric(quartiles[num_pair[2]])
  lower_pair <- c(lower_pair, min(q1, q2))
  upper_pair <- c(upper_pair, max(q1, q2))
}
table(lower_pair)
table(upper_pair)
# negatively correlated stuff tends to be lower abundance + higher abundance

lower_pair <- c()
upper_pair <- c()
for(tax_pair in top_pos) {
  num_pair <- as.numeric(str_split(labels.reordered[tax_pair], "_")[[1]])
  q1 <- as.numeric(quartiles[num_pair[1]])
  q2 <- as.numeric(quartiles[num_pair[2]])
  lower_pair <- c(lower_pair, min(q1, q2))
  upper_pair <- c(upper_pair, max(q1, q2))
}
table(lower_pair)
table(upper_pair)

# positives correlated stuff tends to be ???

