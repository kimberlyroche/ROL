library(ROL)
library(driver)
library(dplyr)
library(ggplot2)
library(stringr)
library(gridExtra)
library(cowplot)

tax_level <- "ASV"
logratio <- "clr"
model_list <- get_fitted_model_list(tax_level=tax_level, MAP=TRUE)
pairs_obj <- get_pairwise_correlations(tax_level=tax_level, logratio="clr")
labels <- pairs_obj$labels
interactions <- pairs_obj$interactions # hosts x pairs

ranked_interactions <- gather_array(interactions, "value", "host", "pair")
ranked_interactions$rank <- as.factor(order(ranked_interactions$value))
ranked_interactions <- ranked_interactions[ranked_interactions$rank,] # min correlation [1] to max [nrow()]
rownames(ranked_interactions) <- 1:nrow(ranked_interactions)

# sample along the rows
negative_idx_bounds <- c(1, max(as.numeric(rownames(ranked_interactions[ranked_interactions$value <= -0.5,]))))
positive_idx_bounds <- c(min(as.numeric(rownames(ranked_interactions[ranked_interactions$value >= 0.5,]))), nrow(ranked_interactions))
idx <- c(round(seq(negative_idx_bounds[1], negative_idx_bounds[2], length.out=6)),
         round(seq(positive_idx_bounds[1], positive_idx_bounds[2], length.out=6)))

# choose indices by quantile instead; there's probably a more elegant way than this
# percentile <- ecdf(ranked_interactions$value)
# ranked_interactions$percentile <- round(sapply(ranked_interactions$value, percentile),2)
# pp <- c(0.0)
# idx <- c()
# for(p in pp) {
#   hits <- ranked_interactions[ranked_interactions$percentile == p,]
#   idx <- c(idx, as.numeric(rownames(hits))[1])
# }

# get taxonomy for readable names
taxonomy <- assign_concise_taxonomy(tax_level=tax_level, logratio="none")

for(ii in 1:length(idx)) {
  i <- idx[ii]
  host <- model_list$hosts[ranked_interactions[i,]$host]
  pair <- pairs_obj$labels[ranked_interactions[i,]$pair]
  value <- ranked_interactions[i,]$value
  microbe_pair <- as.numeric(strsplit(pair, "_")[[1]])
  # assumes tax level is ASV, logratio is CLR
  microbe_pair_labels <- c(paste0("CLR(ASV in ",taxonomy[microbe_pair[1]],")"),
                           paste0("CLR(ASV in ",taxonomy[microbe_pair[2]],")"))
  cat(paste0("Evaluating pair (",microbe_pair[1],") ",microbe_pair_labels[1]," vs. (",microbe_pair[2],") ",
    microbe_pair_labels[2]," in host ",host,"...\n"))

  fit_filename <- file.path("output","model_fits",tax_level,paste0(host,"_bassetfit.rds"))
  fit_obj <- read_file(fit_filename)

  # a dumb hack for now; these need to be global
  SE_sigma <<- fit_obj$kernelparams$SE_sigma
  SE_rho <<- fit_obj$kernelparams$SE_rho
  PER_sigma <<- fit_obj$kernelparams$PER_sigma

  # note: these predictions are in the ALR
  predict_obj <- make_posterior_predictions(fit_obj$X, fit_obj$fit)

  # convert representation if necessary
  Eta <- predict_obj$Y_predict
  if(logratio != "alr") {
    lr_ys <- alr(t(fit_obj$Y) + 0.5)
    if(logratio == "clr") {
      Eta <- alrInv_array(Eta, fit_obj$fit$D, 1)
      Eta <- clr_array(Eta, 1)
      lr_ys <- clr(t(fit_obj$Y) + 0.5)
    }
    if(logratio == "ilr") {
      Eta <- alrInv_array(Eta, fit_obj$fit$D, 1)
      Eta <- ilr_array(Eta, 1)
      lr_ys <- ilr(t(fit_obj$Y) + 0.5)
    }
  } 
  
  # get mean trajectory for microbe 1
  posterior_samples.1 <- gather_array(Eta[microbe_pair[1],,], "logratio_value", "observation", "sample_number")
  post_quantiles.1 <- as.data.frame(posterior_samples.1 %>%
    group_by(observation) %>%
    summarise(mean = mean(logratio_value)) %>%
    ungroup())
  # center; we'll plot on top of each other (we're just interested in how/whether "dynamics" line up)
  post_quantiles.1$mean <- scale(post_quantiles.1$mean, scale=F)

  posterior_samples.2 <- gather_array(Eta[microbe_pair[2],,], "logratio_value", "observation", "sample_number")
  post_quantiles.2 <- as.data.frame(posterior_samples.2 %>%
    group_by(observation) %>%
    summarise(mean = mean(logratio_value)) %>%
    ungroup())
  post_quantiles.2$mean <- scale(post_quantiles.2$mean, scale=F)

  df <- data.frame(x=c(post_quantiles.1$observation, post_quantiles.2$observation),
                   y=c(post_quantiles.1$mean, post_quantiles.2$mean),
                   label=c(rep(paste0(microbe_pair_labels[1]," (#",microbe_pair[1],")"), nrow(post_quantiles.1)),
                           rep(paste0(microbe_pair_labels[2]," (#",microbe_pair[2],")"), nrow(post_quantiles.2))))

  p.1 <- ggplot(df) +
    geom_line(aes(x=x, y=y, color=label)) +
    xlab("days from first sample") +
    ylab("centered mean CLR abundance") +
    ggtitle(paste0("Correlation: ",round(value, 2))) +
    theme(legend.position="bottom")
  if(sign(value) < 0) {
    append_str = "neg"
  } else {
    append_str = "pos"
  }

  df2 <- data.frame(x=post_quantiles.1$mean, y=post_quantiles.2$mean)

  p.2 <- ggplot(df2) +
    geom_point(aes(x=x, y=y)) +
    xlab("centered mean CLR abundance 1") +
    ylab("centered mean CLR abundance 2")

  p <- plot_grid(p.1, p.2, align="h", nrow=1, rel_widths=c(3/4, 1/4))
  ggsave(paste0("output/sanity_check_pair_",str_pad(ii, 2, pad="0"),".png"),
    p, dpi=100, units="in", height=3, width=12)

}
