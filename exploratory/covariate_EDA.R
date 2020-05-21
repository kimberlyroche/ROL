library(tidyverse)
library(phyloseq)
library(ggplot2)
library(pracma)

data <- readRDS("input/ps_w_covs.rds")
metadata <- sample_data(data)
# names(metadata)

# baboons have numbers (e.g. "Baboon_562" instead of short names in this data set)

# (1) select individuals with at least 40 samples
over40_hosts <- as.character(unique(pull(metadata %>%
  select(c("host", "collection_date")) %>%
  group_by(host) %>%
  mutate(count = n()) %>%
  filter(count >= 40) %>%
  select(host), host)))

md.partial <- subset_samples(metadata, host %in% over40_hosts)

if(FALSE) {
  # (2) quick visualizations
  
  quick_viz <- function(md, feature, legend = FALSE) {
    p <- ggplot(md) +
      geom_point(aes_string(x = "diet_PC1", y = "diet_PC2", color = feature))
    if(feature != "host" & feature != "season") {
      # i.e. if feature is continuous
      p <- p + scale_color_gradient(low = "gray", high = "red")
    }
    if(!legend) {
      p <- p + theme(legend.position = "none")
    }
    suppressWarnings(show(p))
  }
  
  temp <- md.partial[sample(1:nrow(md.partial))[1:500],]
  
  quick_viz(temp, "host")
  quick_viz(temp, "season", legend = TRUE)
  
  # climate covariates we might consider: hydro_year, season, rain_monthly, rain_annual, tempmax_monthly, tempmax_annual
  quick_viz(temp, "rain_monthly", legend = TRUE) # not a lot of variation in this feature; rainfall mostly low
  quick_viz(temp, "tempmax_monthly", legend = TRUE)
  
  # other interesting stuff: no obvious relationship between location and diet
  quick_viz(temp, "lat", legend = TRUE)
  quick_viz(temp, "lat", legend = TRUE)

  # (2) do individuals "travel" around the space of diet a lot?
  prop_travel <- function(md, feature) {
    min.feat <- min(md.partial[[feature]], na.rm = TRUE)
    max.feat <- max(md.partial[[feature]], na.rm = TRUE)
    span.feat <- abs(max.feat - min.feat)
    var_sym <- sym(feature)
    host_extrema <- as.data.frame(suppressWarnings(md.partial %>%
      select(c("host", feature)) %>%
      group_by(host) %>%
      summarize(min_PC = min(!!var_sym, na.rm = TRUE), max_PC = max(!!var_sym, na.rm = TRUE))))
    host_extrema <- data.matrix(host_extrema[,2:3])
    # show proportion of travel on diet_PC1 for hosts
    hist(apply(host_extrema, 1, function(x) (abs(x[1] - x[2])/span.feat)), xlim=c(0,1), main=feature)
  }
  
  # what's with the NA's in diet?
  prop_travel(md.partial, "diet_PC1")
  prop_travel(md.partial, "diet_PC2")
  prop_travel(md.partial, "rain_monthly")
  prop_travel(md.partial, "tempmax_monthly")
}

# (3) autocorrelation; what's the scale and decay of diet (etc.)
viz_ac <- function(md, feature, host_limit  = NULL) {
  hosts <<- unique(md)$host
  ac <- data.frame(host = c(), diff = c(), val1 = c(), val2 = c())
  if(is.null(host_limit)) {
    host_limit <- length(hosts)
  }
  use_hosts <- hosts[sample(1:length(hosts))[1:host_limit]]
  for(host in use_hosts) {
    host <<- host
    cat("Host:",host,"\n")
    feature <<- feature
    temp <- subset_samples(md, host == hosts[1])
    temp <- temp[,c(feature, "collection_date")]
    # collection dates look already to be in order
    for(i in 1:(nrow(temp)-1)) {
      for(j in (i+1):nrow(temp)) {
        # only evaluate up to 2 years; data is pretty thin after than
        time.i <- temp$collection_date[i]
        time.j <- temp$collection_date[j]
        diff.ij <- round(difftime(time.j, time.i, units = "weeks"))
        if(diff.ij < 105) {
          val.i <- temp[[feature]][i]
          val.j <- temp[[feature]][j]
          if(!is.na(val.i) & !is.na(val.j)) {
            ac <- rbind(ac, data.frame(host = as.character(host), diff = as.numeric(diff.ij), val1 = val.i, val2 = val.j))
          }
        }
      }
    }
  }
  
  # exclude differences with no variance in the measurements
  # this appears to happen for 1 or 2 of the long lags where every instance of the lag across individuals comes
  # from samples in the same pair of months (e.g. 2001-08 diet x 2008-03 diet)
  viable_differences <- pull(ac %>%
    group_by(diff) %>%
    summarise(std1 = sd(val1), std2 = sd(val2)) %>%
    filter(std1 > 0 & std2 > 0) %>%
    select(diff), "diff")
  
  # calculation correlation at each difference
  temp.2 <- ac[ac$diff %in% viable_differences,] %>%
    group_by(diff) %>%
    summarise(ac = cor(val1, val2)) %>%
    select(diff, ac)
  
  df <- as.data.frame(temp.2)
  p <- ggplot(df) +
    geom_line(aes(x = diff, y = ac)) +
    geom_point(aes(x = diff, y = ac)) +
    xlab("lag (weeks)") +
    ylab("correlation")
  ggsave(paste0(feature,"_autocorrelation.png"), p, units = "in", dpi = 100, height = 6, width = 10)
}

host_limit = 50
#viz_ac(md.partial, feature = "diet_PC1", host_limit = host_limit)
#viz_ac(md.partial, feature = "rain_monthly", host_limit = host_limit)
viz_ac(md.partial, feature = "tempmax_monthly", host_limit = host_limit)

if(FALSE) {
  PER <- function(X, sigma=1, rho=1, period=24, jitter=0) {
    dist <- as.matrix(dist(t(X)))
    G <- sigma^2 * exp(-2*(sin(pi*dist/period)^2)/(rho^2)) + jitter*diag(ncol(dist))
    return(G)
  }
  
  # N <- 365
  # X <- matrix(1:N, 1, N)
  # G <- PER(X, sigma = 1, rho = 1, period = 365) # larger rho is wider bandwidth, longer lasting correlation
  # plot(G[1,])
  
  # build a periodic kernel for diet
  
  # diet kernel; diet PCs are columns 24 - 36
  md.partial.2 <- subset_samples(md.partial, hosts %in% c("Baboon_1", "Baboon_2"))
  X <- md.partial.2[,24:26]
  kern <- PER(t(X), sigma = 1, rho = 1, period = abs(max(X) - min(X)))
  image(kern)
  
  # rainfall
  X <- md.partial.2[,"rain_monthly",drop=F]
  kern <- PER(t(X), sigma = 1, rho = 0.5, period = abs(max(X) - min(X)))
  rownames(kern) <- NULL
  colnames(kern) <- NULL
  round(kern[1:15,1:15], 2)
  image(kern)
  
  # tempmax
  X <- md.partial.2[,"tempmax_monthly",drop=F]
  kern <- PER(t(X), sigma = 1, rho = 1, period = abs(max(X) - min(X)))
  rownames(kern) <- NULL
  colnames(kern) <- NULL
  round(kern[1:15,1:15], 2)
  image(kern)
}

# try fitting DUIKER's samples with and without these extra kernels

# "Baboon_103" must be DUIKER at 181 samples

library(ROL)
data <- load_data(tax_level="ASV")
tax_level <- "ASV"
n_samples <- 100
MAP <- FALSE
alr_ref <- NULL

# global assign is a hack seemingly necessary for this phyloseq::subset_samples function call
host <<- "DUI"
host_data <- subset_samples(data, sname==host)

# append these covariates
md.DUI <- md.partial[md.partial$host == "Baboon_103",]

# encode observations as differences from baseline in units of days
host_metadata <- sample_data(host_data)
baseline_date <- host_metadata$collection_date[1]
observations <- sapply(host_metadata$collection_date, function(x) round(difftime(x, baseline_date, units="days"))) + 1

# pull out the count table
Y <- t(otu_table(host_data)@.Data)
dim(observations) <- c(1, length(observations))

# get dimensions
D <- nrow(Y)
N <- ncol(Y)
  
# stray uses the D^th element as the ALR reference by default
# if we'd like to use a different reference, do some row shuffling in Y to put the reference at the end
if(!is.null(alr_ref)) {
  Y <- Y[c(setdiff(1:D,alr_ref),alr_ref),]
}

# strip off sequence variant labels
colnames(Y) <- NULL
rownames(Y) <- NULL

# back-calculate the squared exponential bandwidth parameter by finding a bandwidth that gives
# a desired minimum correlation at the number of days specified by SE_days_to_baseline
SE_sigma <- 2.2
PER_sigma <- 0.5
dc <- 0.1 # desired minimum correlation
SE_days_to_baseline <- 90
SE_rho <- sqrt(-SE_days_to_baseline^2/(2*log(dc))) # back calculate the decay

# define the composite kernel over samples
Gamma.extra <- function(X) {
  jitter <- 1e-08
  SE(X[1,,drop=F], sigma=SE_sigma, rho=SE_rho, jitter=jitter) +
    PER(X[1,,drop=F], sigma=(PER_sigma/2), rho=1, period=365, jitter=jitter) +
    PER(X[2:4,,drop=F], sigma=(PER_sigma/2), rho=1, period=abs(max(X[2:4,])-min(X[2:4])), jitter=jitter) #+
#    diag(ncol(X))*0.1
}

X.extra <- observations
X.extra <- rbind(X.extra, t(md.DUI[,24:26]))
#X.extra <- rbind(X, t(md.DUI[,"rain_monthly",drop=F]))
#X.extra <- rbind(X, t(md.DUI[,"tempmax_monthly",drop=F]))

# use a default ALR prior
prior_params <- default_ALR_prior(D)
  
# define the prior over baselines
alr_ys <- driver::alr((t(Y) + 0.5))
alr_means <- colMeans(alr_ys)
Theta <- function(X) matrix(alr_means, D-1, ncol(X))
  
# MAP
n_samples <- 10

fit.extra <- stray::basset(Y, X.extra, prior_params$upsilon, Theta, Gamma.extra, prior_params$Xi,
                           n_samples=n_samples, ret_mean=MAP)


# define the composite kernel over samples
Gamma.basic <- function(X) {
  SE(X[1,,drop=F], sigma=SE_sigma, rho=SE_rho, jitter=1e-08) +
    PER(X[1,,drop=F], sigma=PER_sigma, rho=1, period=365, jitter=1e-08)
}

X.basic <- observations

fit.basic <- stray::basset(Y, X.basic, prior_params$upsilon, Theta, Gamma.basic, prior_params$Xi,
                           n_samples=n_samples, ret_mean=MAP)

fit.extra$logMarginalLikelihood # a tiny bit better
fit.basic$logMarginalLikelihood

crude_mean_Sigma.extra <- apply(fit.extra$Sigma, c(1,2), mean)
crude_mean_Sigma.basic <- apply(fit.basic$Sigma, c(1,2), mean)

image(crude_mean_Sigma.extra)
image(crude_mean_Sigma.basic)

# visualize posterior intervals over Lambda, Eta
posterior_plot <- function(X_predict, Y_predict, taxon_idx) {
  sample_set <- Y_predict[taxon_idx,,]
  sample_df <- gather_array(sample_set, "logratio_value", "observation", "sample_number")
  quantiles <- sample_df %>%
    group_by(observation) %>%
    summarise(p2.5 = quantile(logratio_value, prob=0.025),
              p25 = quantile(logratio_value, prob=0.25),
              mean = mean(logratio_value),
              p75 = quantile(logratio_value, prob=0.75),
              p97.5 = quantile(logratio_value, prob=0.975)) %>%
    ungroup()
  
  lr_tidy <- data.frame(observation = c(observations), logratio_value = alr_ys[,taxon_idx])

  p <- ggplot(quantiles, aes(x=observation, y=mean)) +
    geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
    geom_ribbon(aes(ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
    geom_line(color="blue") +
    geom_point(data=lr_tidy, aes(x=observation, y=logratio_value), alpha=0.5) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle=45)) +
    ylab("LR coord")
  show(p)
}

X_predict.basic <- t(1:(max(X.basic)))
predicted.basic <- predict(fit.basic, X_predict.basic, response="Eta", iter=fit.basic$iter)
posterior_plot(X_predict.basic, predicted.basic, 1)

# interpolate
# n <- 100
# x <- 1:n
# y <- 0.5*x + rnorm(n, 0, 1)
# plot(x,y)
# censor <- as.logical(rbinom(n, 1, 0.5))
# x.censor <- x[!censor]
# y.censor <- y[!censor]
# plot(x.censor, y.censor)
# pair <- approx(x.censor, y.censor, xout = 1:n) # warning
# lines(pair$x, pair$y, type="p", col="red")

# need to interpolate diet PCs
n <- max(X.extra[1,,drop=F])
# X_predict.extra <- t(1:n)
nn <- 700
X_predict.extra <- t(1:nn)
trunc_obs <- c(observations)
select_idx <- trunc_obs <= nn
trunc_obs <- trunc_obs[select_idx]
# X_predict.extra <- rbind(X_predict.extra, pchip(xi = trunc_obs, yi = X.extra[2,select_idx], x = 1:nn))
# X_predict.extra <- rbind(X_predict.extra, pchip(xi = trunc_obs, yi = X.extra[3,select_idx], x = 1:nn))
# X_predict.extra <- rbind(X_predict.extra, pchip(xi = trunc_obs, yi = X.extra[4,select_idx], x = 1:nn))
# interpolate the mean
yi.2 <- rep(mean(X.extra[2,select_idx]), ncol(X_predict.extra))
yi.2[trunc_obs] <- X.extra[2,select_idx]
X_predict.extra <- rbind(X_predict.extra, yi = yi.2)
yi.3 <- rep(mean(X.extra[3,select_idx]), ncol(X_predict.extra))
yi.3[trunc_obs] <- X.extra[3,select_idx]
X_predict.extra <- rbind(X_predict.extra, yi = yi.3)
yi.4 <- rep(mean(X.extra[4,select_idx]), ncol(X_predict.extra))
yi.4[trunc_obs] <- X.extra[4,select_idx]
X_predict.extra <- rbind(X_predict.extra, yi = yi.4)

plot(X_predict.extra[2,])

# kernel is not PSD over big prediction range
# TO DO: fiddle with sigma and rho to fix this

test <- PER(X_predict.extra[2,,drop=F], sigma=1, rho=1, period=abs(max(X_predict.extra[2:4,])-min(X_predict.extra[2:4])), jitter=1e-08)
vals <- eigen(test)$values
min(vals)

plot(test[1,])

fit.extra$Gamma <- Gamma.extra
test <- fit.extra$Gamma(X_predict.extra)
vals <- eigen(test)$values
min(vals)
plot(test[1,])

image(test)

predicted.extra <- predict(fit.extra, X_predict.extra, response="Eta", iter=fit.extra$iter)

posterior_plot(X_predict.basic, predicted.basic, 1)




# if(MAP) {
#   fit <- fix_MAP_dims(fit)
# }

# pack up results
# fit_obj <- list(Y=Y, alr_ys=alr_ys, alr_ref=alr_ref, X=observations, fit=fit)
  
# a hack: for later prediction, these will need to be available in the workspace for Gamma
# fit_obj$kernelparams$SE_sigma <- SE_sigma
# fit_obj$kernelparams$SE_rho <- SE_rho
# fit_obj$kernelparams$PER_sigma <- PER_sigma






