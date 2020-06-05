library(stray)

# load the stuff that needs to be global 
data <- load_data(tax_level="ASV")
params <- formalize_parameters(data)
data <- subset_samples(data, sname == host)
metadata <- sample_data(data)

# read diet covariate data
data.diet <- readRDS("input/ps_w_covs.RDS")
data.name_mapping <- read.csv("input/host_subject_id_to_sname_key.csv")
data.name_mapping <- unique(data.name_mapping[,c("sname","host_subject_id2")])
host.num <- as.character(data.name_mapping[data.name_mapping$sname == host,]$host_subject_id2)
data.diet <- subset_samples(data.diet, host == host.num)
metadata.diet <- sample_data(data.diet)

# pull out the count table (taxa x samples)
Y <- t(otu_table(data)@.Data)

# pull dimensions (again)
D <- nrow(Y)
N <- ncol(Y)

colnames(Y) <- NULL
rownames(Y) <- NULL

Y <- Y[c(setdiff(1:D,params$alr_ref),params$alr_ref),]

# define the composite kernel over samples
Gamma.basic <- function(X) {
  dc <- 0.1 # desired minimum correlation
  rho <- sqrt(-90^2/(2*log(dc))) # back calculate the decay (90 days to a drop-off of)
  Gamma_scale <- params$total_variance
  SE(X, sigma = sqrt(Gamma_scale), rho = rho, jitter = 1e-08)
}

# define the composite kernel over samples
# FYI (myself): the optimization is pretty brittle here to the scale of Gamma giving failed (?) convergence
#               (Cholesky of Hessian fails)
Gamma.extra <- function(X) {
  jitter <- 1e-08
  # back calculate the decay in correlation to approx. 0.1 at 90 days
  dc <- 0.1 # desired minimum correlation
  rho <- sqrt(-90^2/(2*log(dc))) # back calculate the decay (90 days to a drop-off of)
  Gamma_scale <- params$total_variance
  base_sigma_sq <- Gamma_scale * 0.5
  PER_sigma_sq <- Gamma_scale * 0.5
  SE(X[1,,drop=F], sigma = sqrt(base_sigma_sq), rho = rho, jitter = jitter) +
    PER(X[1,,drop=F], sigma = sqrt(PER_sigma_sq/2), rho = 1, period = 365, jitter = jitter) +
    SE(X[2:6,,drop=F], sigma = sqrt(PER_sigma_sq/2), rho = 1, jitter = jitter) # a few diet PCs
}

# ALR prior for Sigma (bacterial covariance)
upsilon <- D - 1 + 10 # specify low certainty/concentration
GG <- cbind(diag(D-1), -1) # log contrast for ALR with last taxon as reference
Xi <- GG%*%(diag(D))%*%t(GG) # take diag as covariance over log abundances
Xi <- Xi*(upsilon-D-1)

# define the prior over baselines as the empirical mean alr(Y)
alr_ys <- driver::alr((t(Y) + 0.5))
alr_means <- colMeans(alr_ys)
Theta <- function(X) matrix(alr_means, D-1, ncol(X))






host <- "DUI"
data <- readRDS(paste0(host,"_fit_basic.rds"))
Lambda_hat <- apply(data$Lambda, c(1,2), mean)
Sigma_hat <- apply(data$Sigma, c(1,2), mean)
Theta_applied <- data$Theta(data$X)
Gamma_applied <- data$Gamma(data$X)
L_Gamma <- t(chol(Gamma_applied))
