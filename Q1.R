library(ROL)
library(matrixsampling)
library(driver)
library(stray)
library(ggplot)
library(RColorBrewer)
library(vegan)
library(dplyr)
library(phyloseq)

# (A) is it correct to validate by drawing from the GP model?
# (B) is it correct to think about \eta \sim N(\Phi, \Sigma, I_N) where \Phi is some specified baseline?

# load an example fitted model
model_list <- get_fitted_model_list(tax_level="ASV", MAP=TRUE)
n_hosts <- length(model_list$hosts)
n_samples <- 10
D <- NULL
hosts <- model_list$hosts[1:n_hosts]
simulated_samples <- list()
for(i in 1:n_hosts) {
  data <- readRDS(model_list$model_list[i])
  host <- hosts[i]
  fitted_model <- data$fit
  if(is.null(D)) {
    D <- fitted_model$D
  }
  fitted_model.clr <- to_clr(fitted_model)
  prop_ys <- alrInv_array(data$alr_ys, coord=2) # returns samples x taxa orientation
  clr_ys <- clr_array(prop_ys, 2)
  baseline <- colMeans(clr_ys)
  # this is a ridiculous means of smoothing over numeric issues with these estimates
  # that give some of the dynamics matrices vanishingly small imaginary eigenvalues
  # or tiny fluctuations that make them not-perfectly-symmetric
  dynamics_sq <- chol(fitted_model.clr$Sigma[,,1] + diag(D)*1e-08)
  dynamics <- t(dynamics_sq)%*%dynamics_sq
  simulated_samples[[host]] <- cbind(baseline,
    rmatrixnormal(1,
      matrix(baseline, D, n_samples),
      dynamics,
      diag(n_samples)
    )[,,1]
  )
}

if(FALSE) {
  # barplot to visualize
  for(i in 1:n_hosts) {
    host <- hosts[i]
    cat("Visualizing host:",host,"\n")
    props <- clrInv_array(simulated_samples[[host]][,1:10], 1)
    prop_df <- driver::gather_array(props, "proportion", "taxon", "sample")
    prop_df$taxon <- as.factor(prop_df$taxon)

    coul = brewer.pal(4, "Spectral")
    coul = colorRampPalette(coul)(length(unique(prop_df$taxon)))

    p.alr <- ggplot(prop_df, aes(x=sample, y=proportion, fill=taxon)) +
      geom_bar(position="fill", stat="identity", width=0.65) +
      scale_fill_manual(values=coul) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.text=element_text(size=8)) +
      theme(legend.position="none") +
      theme_light() +
      theme(legend.position = "none")
    ggsave(paste0("output/simulated_samples_",host,".png"), p.alr, dpi=150, units="in", height=3, width=6)
  }
}

sample_labels <- c()
embed_matrix <- matrix(NA, n_hosts*(n_samples + 1), D)
for(i in 1:n_hosts) {
  host <- hosts[i]
  sample_labels <- c(sample_labels, rep(host, n_samples + 1))
  for(j in 1:(n_samples+1)) {
    offset <- (i-1)*(n_samples+1) + j
    embed_matrix[offset,] <- simulated_samples[[host]][,j]
  }
}

distances <- dist(embed_matrix)

# visualize via embedding
embedding <- cmdscale(distances, k=2, eig=TRUE)

df <- data.frame(x=embedding$points[,1], y=embedding$points[,2], host=sample_labels)
p <- ggplot(df) +
  geom_point(aes(x=x, y=y, color=host)) +
  xlab("PCA 1") +
  ylab("PCA 2")
ggsave("output/simulated_samples_embedding.png", p, dpi=150, units="in", height=6, width=9)

# run PERMANOVA
d <- data.frame(labels=sample_labels)
obj <- adonis(distances ~ labels, data=d, permutations=1000)

# report r-squared
R2 <- obj$aov.tab$R2
cat("Host r-squared:",round(R2[1], 3),"\n")

data <- readRDS("input/filtered_ASV_5_20.rds")
metadata <- sample_data(data)
primary_group <- suppressWarnings(metadata %>%
                                 select(c("sname", "collection_date", "grp")) %>%
                                 group_by(sname, grp) %>%
                                 tally() %>%
                                 slice(which.max(n)))

group_labels <- as.character(d$labels)
for(i in 1:length(group_labels)) {
  group_labels[i] <- primary_group[primary_group$sname == group_labels[i],]$grp
}
group_labels <- as.factor(group_labels)

group_labels <- data.frame(labels=group_labels)

# run PERMANOVA; >1000 permutations doesn't seem to change the result
obj <- adonis(distances ~ labels, data=group_labels, permutations=1000)

# report r-squared
R2 <- obj$aov.tab$R2
cat("Group r-squared:",round(R2[1], 3),"\n")

