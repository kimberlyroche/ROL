library(ROL)
library(phyloseq)
library(ggplot2)

# visualizations for Beth's R01

tax_level <- "family"
data <- load_data(tax_level=tax_level, count_threshold=5, sample_threshold=0.2)

# render barplots for individuals who at disparate points in dynamics ordination
hosts <- list(DUN = "Dunlin",
              LEB = "Lebanon",
              HON = "Honey",
              WIP = "Wiper")
for(host in names(hosts)) {
    cat(paste0("Visualizing host ",host,"...\n"))
    plot_timecourse(data, host = host)
}

# render family-level embedding (coord 1 x 2) and label just these hosts' posterior samples

# load family-level full posterior ordination and labels
ordination <- readRDS("output/plots/family/Sigma_ordination.rds")
ordination$labels <- as.character(ordination$labels)
# rename to full names
for(host in names(hosts)) {
    ordination[ordination$labels == host,]$labels <- hosts[[host]]
}
# clear the rest
ordination[!(ordination$labels %in% hosts),]$labels <- NA

df <- data.frame(x = ordination[ordination$coord == 1,]$value,
                 y = ordination[ordination$coord == 2,]$value,
                 Host = as.factor(ordination[ordination$coord == 1,]$labels))


p <- ggplot(df, aes(x = x, y = y, color = Host)) +
    geom_point() +
    scale_color_discrete(na.value="#bbbbbb") +
    guides(fill = guide_legend(title = "Host name")) +
    xlab("PCoA 1") +
    ylab("PCoA 2")
ggsave("output/plots/ordination_R01.png", p, dpi = 100, units = "in", width = 10, height = 8)

