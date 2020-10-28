library(phyloseq)
library(ggplot2)
library(GGally) # for ggpairs

# Some covariates have missing values; we'll imput the mean for this as a quick
# and dirty way of dealing with them
impute_mean <- function(values) {
  values[is.na(values)] <- mean(values, na.rm = TRUE)
  values
}

# Pull data/metadata that includes environmental covariates
data <- readRDS("C:/Users/kim/Documents/ROL/exploratory/input/ps_w_covs.RDS")
md <- sample_data(data)

# Show covariate names
names(md)

# Get hosts names (e.g. "Baboon_###") ordered by max to min samples (down to minimum 75)
well_sampled_hosts <- names(which(sort(table(md$host), decreasing = TRUE) >= 75))

# Pull a host
host_label <- well_sampled_hosts[1]

# Pull interesting covariates
rain_monthly <- md$rain_monthly[md$host == host_label]
tempmax_monthly <- md$tempmax_monthly[md$host == host_label]
diet_PC1 <- md$diet_PC1[md$host == host_label]
diet_PC2 <- md$diet_PC2[md$host == host_label]
diet_PC3 <- md$diet_PC3[md$host == host_label]
diet_PC4 <- md$diet_PC4[md$host == host_label]
diet_PC5 <- md$diet_PC5[md$host == host_label]
diet_PC6 <- md$diet_PC6[md$host == host_label]

# Get rid of NAs the dumb way
rain_monthly <- impute_mean(rain_monthly)
tempmax_monthly <- impute_mean(tempmax_monthly)
diet_PC1 <- impute_mean(diet_PC1)
diet_PC2 <- impute_mean(diet_PC2)
diet_PC3 <- impute_mean(diet_PC3)
diet_PC4 <- impute_mean(diet_PC4)
diet_PC5 <- impute_mean(diet_PC5)
diet_PC6 <- impute_mean(diet_PC6)

df <- data.frame(rain = rain_monthly,
                 temp = tempmax_monthly,
                 diet1 = diet_PC1,
                 diet2 = diet_PC2,
                 diet3 = diet_PC3,
                 diet4 = diet_PC4,
                 diet5 = diet_PC5,
                 diet6 = diet_PC6)
ggpairs(df)
