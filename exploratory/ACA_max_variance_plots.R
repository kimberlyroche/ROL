library(phyloseq)

# What does the real data look like?
data <- load_data(tax_level = "ASV")
ACA_data <- subset_samples(data, sname == "ACA")
ACA_date <- sample_data(ACA_data)$collection_date
ACA_season <- sample_data(ACA_data)$season
baseline_date <- min(ACA_date)
days <- unname(round(sapply(ACA_date, function(x) difftime(x, baseline_date, units = "days"))) + 1)
tax <- tax_table(ACA_data)
counts <- otu_table(ACA_data)@.Data
unname(colMeans(counts))
hist(counts, breaks = 20)
clr_counts <- clr_array(counts + 0.1, parts = 2)
hist(clr_counts[,sample(1:ncol(clr_counts))[1]], breaks = 20)
dim(counts)  
apply(clr_counts, 2, var)

# The upper limit of variation is pretty big
tax_label <- tax@.Data[29,6]
df <- data.frame(clr_counts = clr_counts[,29], day = days, season = as.factor(ACA_season))
ggplot(df, aes(x = day, y = clr_counts, color = season)) +
  geom_point() +
  ylab(paste0("unmodeled CLR abundance genus ",tax_label))

df <- data.frame(counts = counts[,29], day = days, season = as.factor(ACA_season))
ggplot(df, aes(x = day, y = counts, color = season)) +
  geom_point() +
  ylab(paste0("unmodeled abundance genus ",tax_label))

