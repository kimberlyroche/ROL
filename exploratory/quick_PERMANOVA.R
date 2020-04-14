library(vegan)
library(phyloseq)
library(dplyr)

# load distances
d <- readRDS("output/Sigma_distance_ASV.rds")

# get distances into their own data.frame
dmat <- d$distance_mat*0.5

# run PERMANOVA; >1000 permutations doesn't seem to change the result
obj <- adonis(dmat ~ host_labels, data=d, permutations=1000)

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

group_labels <- d$host_labels
for(i in 1:length(group_labels)) {
	group_labels[i] <- primary_group[primary_group$sname == group_labels[i],]$grp
}
group_labels <- as.factor(group_labels)

group_labels <- data.frame(labels=group_labels)

# run PERMANOVA; >1000 permutations doesn't seem to change the result
obj <- adonis(dmat ~ labels, data=group_labels, permutations=1000)

# report r-squared
R2 <- obj$aov.tab$R2
cat("Group r-squared:",round(R2[1], 3),"\n")
