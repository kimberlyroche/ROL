#' Render bar plots for an individual host's series
#' 
#' @param data a phyloseq object
#' @param host host short name (e.g. ACA)
#' @param gapped render true sample spacing (approximately)
#' @param legend render legend
#' @param show_plot show() plot in addition to rendering it to a file
#' @param selected_samples list of sample IDs to highlight (if visualizing a selection of a host's samples)
#' @return NULL
#' @import phyloseq
#' @import dplyr
#' @import RColorBrewer
#' @import ggplot2
#' @export
#' @examples
#' plot_timecourse(data, host="ACA", gapped=FALSE, legend=TRUE, selected_samples=NULL)
plot_timecourse <- function(data, host=NULL, gapped=FALSE, legend=TRUE, selected_samples=NULL, show_plot=FALSE) {
  if(is.null(host)) {
    cat("Error: omitted host short name in plot_timecourse()!\n")
    return(NULL)
  }
  # global assign is a hack seemingly necessary for this phyloseq::subset_samples function call
  host <<- host
  # subset to selected host
  subsetted_data <- subset_samples(data, sname == host)
  # transform to proportions
  proportions <- transform_sample_counts(subsetted_data, function(x) x / sum(x)) # samples x taxa
  N <- phyloseq::nsamples(proportions)
  cat(paste0("Plotting ",N,"-sample timecourse\n"))
  proportions_df <- psmelt(proportions)
  trimmed_df <- bind_cols(list(OTU=proportions_df$OTU,
                               Sample=proportions_df$sample_id,
                               Abundance=proportions_df$Abundance,
                               SID=proportions_df$sid))
  if(!is.null(selected_samples)) {
    trimmed_df$alpha <- 0.25
    trimmed_df[trimmed_df$SID %in% selected_samples,]$alpha <- 1
  }
  
  # for gapped plots, this is the OTU ID placeholder we'll use for dummy data
  na.string <- ".N/A"
  
  # replace sample ID's with their dates for readability
  metadata <- sample_data(subsetted_data)
  for(i in 1:nrow(trimmed_df)) {
    if(gapped) {
      trimmed_df$Sample[i] <- metadata[metadata$sample_id==trimmed_df$Sample[i],"collection_date"][[1]]
    } else {
      trimmed_df$Sample[i] <- paste(metadata[metadata$sample_id==trimmed_df$Sample[i],"collection_date"][[1]],
                                    metadata[metadata$sample_id==trimmed_df$Sample[i],"sid"][[1]])
    }
  }
  
  # make sure the dates are in order and fix the order by converting to factors
  reordered_df <- trimmed_df[order(trimmed_df$Sample),]
  if(gapped) {
    reordered_df <- reordered_df[,c("OTU","Sample","Abundance")] # kill SID for this case
    # insert empty samples where gaps of 2 weeks or more exist
    gap.days <- 13
    dates_present <- unique(reordered_df$Sample)
    for(d in 1:(length(dates_present)-1)) {
      diff <- as.Date(dates_present[d+1]) - as.Date(dates_present[d])
      next.date <- as.Date(dates_present[d])
      attr(diff, "units") <- "days"
      while(diff > gap.days) {
        next.date <- next.date + gap.days
        reordered_df <- rbind(reordered_df, list(OTU=na.string, Sample=as.character(next.date), Abundance=1))
        diff <- as.Date(dates_present[d+1]) - next.date
      }
    }
    reordered_df <- reordered_df[order(reordered_df$Sample),]
  }
  
  # this label appears to need to be a factor to maintain its ordering
  reordered_df$Sample <- factor(reordered_df$Sample, levels=unique(reordered_df$Sample))

  # replace ASV sequences with their (abbrev.) family (=5) level taxonomy for readability
  img_height = 4
  if(legend) {
    for(i in 1:dim(reordered_df)[1]) {
      # show labels as order/family/genus
      # species is NA for all
      if(reordered_df$OTU[i] != na.string) {
        reordered_df$OTU[i] <- paste(as.vector(tax_table(data)[reordered_df$OTU[i],5]),collapse="/")
      }
    }
    img_height = 6
  }
  
  categories <- unique(reordered_df$OTU)
  coul = brewer.pal(4, "Spectral")
  coul = colorRampPalette(coul)(length(unique(reordered_df$OTU)))
  if(gapped | !legend) {
    coul[1] <- "#DDDDDD"
  }
  # img_width <- round(N/4)
  # fix the dimensions of these at 4 x 12 in.
  if(!is.null(selected_samples)) {
    p <- ggplot(reordered_df, aes(x=Sample, y=Abundance, fill=OTU, alpha=alpha))
  } else {
    p <- ggplot(reordered_df, aes(x=Sample, y=Abundance, fill=OTU))
  }
  p <- p + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=coul) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.text=element_text(size=8))
  if(legend) {
    p <- p + theme(legend.position="bottom")
  } else {
    p <- p + theme(legend.position="none")
  }
  save_dir <- check_output_dir(c("output","plots"))
  if(show_plot) {
    show(p)
  }
  ggsave(file.path(save_dir,paste0(host,"_timecourse.png")),
         plot=p, scale=1.5, dpi=100, width=12, height=img_height, units="in")
}

#' Calculate lagged autocorrelation
#' 
#' @param data a phyloseq object
#' @param lag_units units of autocorrelation lag ("weeks", "months", "seasons") 
#' @param lag_max maximum lagged difference to evaluate
#' @param use_lr log ratio representation ("ilr", "alr", "clr")
#' @param alr_ref if using the additive log ratio (ALR), which reference taxon index to use
#' @param resample use resampling for uncertainty estimation
#' @param resample_rate percent of samples to randomly subsample if resampling
#' @return data.frame of (bounded, if resampled) autocorrelation, indexed by lag in lag_units
#' @import ggplot2
#' @import driver
#' @export
#' @examples
#' lagged_ac <- calc_autocorrelation(data, lag_units="weeks", lag_max=52, use_lr="ilr", alr_ref=NULL, resample=FALSE, resample_rate=0.2)
calc_autocorrelation <- function(data, lag_units="weeks", lag_max=52, use_lr="ilr", alr_ref=NULL, resample=FALSE, resample_rate=0.2) {
  count_table <- otu_table(data)
  metadata <- sample_data(data)
  
  # select hosts to loop
  hosts <- unique(metadata$sname)

  if(lag_units == "season") {
    # if using season, manually use these YYYYMM as wet-dry season cutoffs
    season_boundaries <- c(200005, 200010, 200105, 200110, 200205, 200210, 200305, 200310,
                           200405, 200410, 200505, 200510, 200605, 200610, 200705, 200710,
                           200805, 200810, 200905, 200910, 201005, 201010, 201105, 201110,
                           201205, 201210, 201305, 201310)
  }
  
  # resampling rounds (if no resampling, 1)
  rounds <- 1
  if(resample) {
    rounds <- 100
  }
  
  lags <- matrix(0, nrow=lag_max, ncol=rounds)

  if(use_lr == "clr") {
    log_ratios <- driver::clr(count_table + 0.5)
  } else if(use_lr == "alr") {
    if(!is.null(alr_ref)) {
      data <- data[c(setdiff(1:nrow(count_table),alr_ref),alr_ref),]
    }
    log_ratios <- driver::alr(count_table + 0.5)
  } else {
    log_ratios <- driver::ilr(count_table + 0.5)
  }
  # normalize
  log_ratios <- scale(log_ratios, center=TRUE, scale=FALSE)

  for(r in 1:rounds) {
    lag_sums <- numeric(lag_max)
    lag_measured <- numeric(lag_max)
    for(host in hosts) {
      # global assign is a hack seemingly necessary for this phyloseq::subset_samples function call
      host <<- host
      # pull sample IDs associated with this individual
      sample_info <- metadata[metadata$sname == host, c("sample_id", "collection_date")]
      # reorder in case not already ordered
      sample_info <- sample_info[order(sample_info$collection_date),]
      # consider only individuals with at least 2 samples (this should already be true)
      if(length(intersect(rownames(log_ratios), sample_info$sample_id)) > 1) {
        sample_lr <- log_ratios[rownames(log_ratios) %in% sample_info$sample_id,,drop=F]
        host_sample_no <- dim(sample_lr)[1]
        do_sample <- rep(1, host_sample_no)
        if(resample) {
          # randomly censor some of the samples
          do_sample <- rbinom(host_sample_no, 1, resample_rate)
        }
        # replace these sample IDs with collection dates, that's what we really want
        # order should be maintained here
        rownames(sample_lr) <- sample_info[sample_info$sample_id %in% rownames(sample_lr), "collection_date"]$collection_date
        
        # get distances between adjacent timepoints in units of lag_units
        d1 <- as.Date(sample_info$collection_date[1])
        time_diff <- matrix(0, nrow=host_sample_no, ncol=host_sample_no)
        if(host_sample_no > 1) {
          for(d1_idx in 1:(host_sample_no-1)) {
            for(d2_idx in (d1_idx+1):host_sample_no) {
              d1 <- as.Date(sample_info$collection_date[d1_idx])
              d2 <- as.Date(sample_info$collection_date[d2_idx])
              time_diff[d1_idx,d2_idx] <- as.numeric(difftime(d2, d1, units="weeks"))
              if(lag_units == "weeks") {
                time_diff[d1_idx,d2_idx] <- ceiling(time_diff[d1_idx,d2_idx])
              } else if(lag_units == "months") {
                time_diff[d1_idx,d2_idx] <- ceiling(time_diff[d1_idx,d2_idx]/4)
              } else if(lag_units == "seasons") {
                d1_ym <- as.numeric(format(d1, "%Y%m"))
                d2_ym <- as.numeric(format(d2, "%Y%m"))
                # distance is number of season boundaries between the samples, indicated in the vector above
                time_diff[d1_idx,d2_idx] <- sum(season_boundaries < d2_ym) - sum(season_boundaries < d1_ym) + 1
              }
            }
          }
        }
        for(lag in 1:lag_max) {
          # find pairs with this lag
          idx <- which(time_diff==lag, arr.ind=T)
          if(dim(idx)[1] >= 1) {
            for(t in 1:dim(idx)[1]) {
              d1_idx <- idx[t,1]
              d2_idx <- idx[t,2]
              if(do_sample[d1_idx] && do_sample[d2_idx]) {
                # given centering, cosine angle == Pearson's correlation
                lag_sums[lag] <- lag_sums[lag] + cor(sample_lr[d1_idx,], sample_lr[d2_idx,])
                lag_measured[lag] <- lag_measured[lag] + 1
              }
            }
          }
          if(lag_sums[lag] > 0 & lag_measured[lag] > 0) {
            lags[lag,r] <- lag_sums[lag]/lag_measured[lag]
          }
        }
      }
    }
  }
  
  if(resample) {
    # get 90% confidence intervals, summarize
    lower_CI <- as.numeric(lag_max)
    upper_CI <- as.numeric(lag_max)
    avg_ac <- as.numeric(lag_max)
    for(lag in 1:lag_max) {
      ac_measured <- sort(lags[lag,])
      lower_CI[lag] <- ac_measured[round(rounds*0.1)]
      upper_CI[lag] <- ac_measured[round(rounds*0.9)]
      avg_ac[lag] <- mean(ac_measured)
    }
    return(data.frame(lag=seq(1,lag_max),
                      ac=avg_ac,
                      lower=lower_CI,
                      upper=upper_CI,
                      units=lag_units))
  } else {
    return(data.frame(lag=seq(1,lag_max),
                      ac=lags,
                      units=lag_units))
  }
}

#' Renders autocorrelation plot (without uncertainty estimates)
#' 
#' @param lags data.frame of lags and estimated autocorrelation (output of calc_autocorrelation())
#' @param show_plot show() plot in addition to rendering it to a file
#' @return NULL
#' @import ggplot2
#' @export
#' @examples
#' lagged_ac <- calc_autocorrelation(data, lag_units="weeks", lag_max=52, use_lr="ilr", alr_ref=NULL, resample=FALSE, resample_rate=0.2)
#' plot_mean_autocorrelation(lagged_ac)
plot_mean_autocorrelation <- function(lagged_ac, show_plot=FALSE) {
  lag_max <- max(lagged_ac$lag)
  p <- ggplot(lagged_ac, aes(x=lag, y=ac)) +
    geom_line(linetype = "dashed") +
    scale_x_continuous(breaks=seq(0,lag_max,1)) +
    scale_y_continuous(breaks=seq(-0.5,1,0.1)) +
    xlab("lag") +
    ylab("ACF") +
    theme_minimal() +
    theme(axis.line.x=element_line(), axis.line.y=element_line()) +
    geom_hline(yintercept = 0)
  save_dir <- check_output_dir(c("output","plots"))
  lag_units <- lagged_ac[1,"units"] # there should be only one unique value
  if(show_plot) {
    show(p)
  }
  ggsave(file.path(save_dir,paste0("autocorrelation_",lag_units,"_",lag_max,".png")),
         plot=p, dpi=100, scale=1.5, width=8, height=4, units="in")
}

#' Renders autocorrelation plot with uncertainty estimates
#' 
#' @param lags data.frame of lags and estimated autocorrelation (output of calc_autocorrelation())
#' @param show_plot show() plot in addition to rendering it to a file
#' @return NULL
#' @import ggplot2
#' @export
#' @examples
#' lagged_ac <- calc_autocorrelation(data, lag_units="weeks", lag_max=52, use_lr="ilr", alr_ref=NULL, resample=TRUE, resample_rate=0.2)
#' plot_bounded_autocorrelation(lagged_ac)
plot_bounded_autocorrelation <- function(lagged_ac, show_plot=FALSE) {
  lag_max <- max(lagged_ac$lag)
  p <- ggplot(lagged_ac, aes(x=lag, y=ac)) +
  geom_line(linetype="solid") +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    scale_x_continuous(breaks=seq(0,lag_max,1)) +
    scale_y_continuous(breaks=seq(-0.5,1,0.1)) +
    xlab("lag") +
    ylab("ACF") +
    theme_minimal() +
    theme(axis.line.x=element_line(), axis.line.y=element_line()) +
    geom_hline(yintercept = 0)
  save_dir <- check_output_dir(c("output","plots"))
  lag_units <- lagged_ac[1,"units"] # there should be only one unique value
  if(show_plot) {
    show(p)
  }
  ggsave(file.path(save_dir,paste0("autocorrelation_",lag_units,"_",lag_max,"_resampled.png")),
         plot=p, dpi=100, scale=1.5, width=8, height=4, units="in")
}


