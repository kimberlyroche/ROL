library(ROL)
library(driver)
library(dplyr)
library(phyloseq)
library(dendextend)

tax_level <- "ASV"

heatmap <- plot_interaction_heatmap(tax_level=tax_level, logratio="clr", Sigmas=NULL,
                                         taxon_idx=NULL, show_plot=FALSE, return_matrix=TRUE)

model_list <- get_fitted_model_list(tax_level=tax_level, MAP=TRUE)

# get host to social group mapping
hosts <<- model_list$hosts
data <- load_data(tax_level)
data <- subset_samples(data, sname %in% hosts)
metadata <- sample_data(data)

primary_group <- suppressWarnings(metadata %>%
                                    select(c("sname", "collection_date", "grp")) %>%
                                    filter(sname %in% hosts) %>% 
                                    group_by(sname, grp) %>%
                                    tally() %>%
                                    slice(which.max(n)))

# create a list indexed by host name
labels <- character(length(hosts))
names(labels) <- hosts
for(host in hosts) {
  labels[host] <- primary_group[primary_group$sname == host,]$grp[[1]]
}
labels <- as.factor(labels)

# cluster and visualize hosts
dd <- dist(heatmap, method = "euclidean")
dmat <- as.matrix(dd)

d_labels <- sapply(hosts, function(x) labels[[x]] )

rownames(dmat) <- d_labels
colnames(dmat) <- d_labels
hc <- hclust(as.dist(dmat))
dend <- as.dendrogram(hc)
labels_colors(dend) <- as.numeric(as.factor(hc$labels))
dend <- set(dend, "labels_cex", 1)

png("output/plots/host_interaction_dendrogram.png", height=1000, width=1000)
plot(dend)
dev.off()

