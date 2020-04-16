library(ROL)
library(driver)
library(ggplot2)
library(phyloseq)
library(dplyr)
library(vegan)

plot_heatmap <- function(mat, filename) {
  df <- gather_array(mat, "sign", "host", "interaction")
  p <- ggplot(df, aes(interaction, host)) +
        geom_tile(aes(fill = sign)) +
        scale_fill_gradient2(low = "darkblue", high = "darkred")
  ggsave(filename, p, units="in", dpi=100, height=5, width=15)
}

tax_level <- "ASV"

calc_posterior_distances(tax_level=tax_level, which_measure="Sigma", which_distance="Riemannian", MAP=TRUE, spike_in=TRUE)

embed_posteriors(tax_level=tax_level, which_measure="Sigma", MAP=TRUE, spike_in=TRUE)

# we'll do our own visualization of the ordination because we're going to switch some things up
coordinates <- read_file(file.path("output","plots",paste0(tax_level,"_MAP"),"Sigma_ordination_spikein.rds"))
centroids <- read_file(file.path("output","plots",paste0(tax_level,"_MAP"),"Sigma_ordination_centroids_spikein.rds"))

hosts <<- setdiff(unique(as.character(centroids$labels)), c("random"))
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
labels[["random"]] <- "0.0"
labels <- as.factor(labels)

df_labels <- sapply(coordinates[coordinates$coord == axis1,]$labels, function(x) labels[[x]])

axis1 <- 1
axis2 <- 2
point_df <- data.frame(ax1=coordinates[coordinates$coord == axis1,]$value,
                       ax2=coordinates[coordinates$coord == axis2,]$value,
                       labels=df_labels)

p <- ggplot() +
  geom_point(data=point_df, aes(x=ax1, y=ax2, color=labels), size=1)
p <- p + theme(legend.position='none')
  
save_dir <- check_output_dir(c("output","plots",paste0(tax_level,"_MAP")))
plot_save_name <- paste0("Sigma_ordination_",axis1,"x",axis2,"_spikein.png")
  
ggsave(file.path(save_dir, plot_save_name), plot=p, dpi=100, width=10, height=6, units="in")

# PERMANOVA
d <- readRDS(paste0("output/Sigma_distance_",tax_level,"_MAP_spikein.rds"))
dmat <- d$distance_mat

# render distance matrix
plot_heatmap(dmat, "spikein_distance_mat.png")

d$host_labels <- sapply(d$host_labels, function(x) if(x != "random") { "true" } else { x } )

# run PERMANOVA; >1000 permutations doesn't seem to change the result
obj <- adonis(dmat ~ host_labels, data=d, permutations=1000)

# report r-squared
R2 <- obj$aov.tab$R2
cat("Host r-squared:",round(R2[1], 3),"\n")

