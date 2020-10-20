library(fido)
library(ROL)
library(phyloseq)

host <- "ACA"
data <- readRDS("input/filtered_family_5_20.rds")
data <- subset_samples(data, sname == host)
md <- sample_data(data)

# pull diet and climate data
data.diet <- readRDS("input/ps_w_covs.RDS")
data.name_mapping <- read.csv("input/host_subject_id_to_sname_key.csv")
data.name_mapping <- unique(data.name_mapping[,c("sname","host_subject_id2")])
host.num <- as.character(data.name_mapping[data.name_mapping$sname == host,]$host_subject_id2)
data.diet <- subset_samples(data.diet, host == host.num)
metadata.diet <- sample_data(data.diet)

# days, weeks, months seem too coarse to be explanatory
# quarters starts to get somewhere but season is still best in the linear model
# days <- md$collection_date
# day0 <- min(days)
# time_baseline <- unname(sapply(days, function(x) difftime(x, day0, units = "weeks")))
# time_baseline <- round(time_baseline / 13) %% 4 # pseudo-halves
season <- as.factor(md$season)
# order is verifiably the same here
diet_PC1 <- scale(metadata.diet$PC1)
rain_monthly <- scale(metadata.diet$rain_monthly)

Y <- t(as.matrix(otu_table(data)@.Data))
rownames(Y) <- NULL
colnames(Y) <- NULL
X <- t(model.matrix(~ season + diet_PC1 + rain_monthly))
D <- nrow(Y)
N <- ncol(Y)
alr_ys <- driver::alr((t(Y) + 0.5))
alr_means <- colMeans(alr_ys)
Theta <- matrix(alr_means, D-1, nrow(X))
var_scale <- 1
upsilon <- D-1+10 # specify low certainty/concentration
GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference
Xi <- GG%*%(diag(D)*var_scale)%*%t(GG) # take diag as covariance over log abundances
Xi <- Xi*(upsilon-D-1)
Gamma <- diag(nrow(X))*var_scale

fit <- pibble(Y, X, Theta = Theta, upsilon = upsilon, Xi = Xi, Gamma = Gamma)
fit.clr <- to_clr(fit)

posterior_summary <- summary(fit.clr, pars="Lambda")$Lambda
focus <- posterior_summary[sign(posterior_summary$p2.5) == sign(posterior_summary$p97.5),]
focus <- unique(focus$coord)
plot(fit.clr, par="Lambda", focus.coord = focus, focus.cov = rownames(X)[2:4])

# Lambda (regression coeffs) seem to capture seasonally varying vs. invariant families
mean_Eta <- apply(fit.clr$Eta, c(1,2), mean)
plot(as.factor(season), mean_Eta[1,]) # seasonal effect
plot(as.factor(season), mean_Eta[15,]) # seasonal effect
plot(as.factor(season), mean_Eta[26,]) # no seasonal effect

# acf(fit.clr$Eta[sample(1:(D-1))[1],,1])
