library(tidyverse)
library(phyloseq)
library(ggplot2)

data <- readRDS("input/ps_w_covs.RDS")
metadata <- sample_data(data)
names(metadata)

# baboons have numbers (e.g. "Baboon_562" instead of short names in this data set)

# (1) select individuals with at least 40 samples
over40_hosts <- as.character(unique(pull(metadata %>%
  select(c("host", "collection_date")) %>%
  group_by(host) %>%
  mutate(count = n()) %>%
  filter(count >= 40) %>%
  select(host), host)))

# "Baboon_103" must be DUIKER at 181 samples

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
  for(host in hosts[1:host_limit]) {
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

host_limit = NULL
viz_ac(md.partial, feature = "diet_PC1", host_limit = host_limit)
#viz_ac(md.partial, feature = "rain_monthly", host_limit = host_limit)
#viz_ac(md.partial, feature = "tempmax_monthly", host_limit = host_limit)

if(FALSE) {
  PER <- function(X, sigma=1, rho=1, period=24, jitter=0){
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








