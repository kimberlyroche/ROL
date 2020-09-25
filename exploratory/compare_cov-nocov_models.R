# Compare model fits with and without covariates in terms of (1) "rug" plots and (2) universality of associations
# Assumes no-covariate MAP fits are in the directory: output/model_fits_nocovariates/ASV_MAP
# Assumes covariate-inclusive MAP fits are in the directory: output/model_fits/ASV_MAP

library(stray)
library(ROL)
library(phyloseq)
library(dplyr)
library(driver)
library(ggplot2)

# Code pulled from universality_score_heatmap.R
calc_xy <- function(vSigmas) {
  if(is.null(dim(vSigmas))) {
    dim(vSigmas) <- c(length(vSigmas), 1)
  }
  shared_direction <- mean(apply(sign(vSigmas), 2, function(z) {
    max(table(z)/length(z))
  }))
  mean_strength <- mean(abs(apply(vSigmas, 2, mean)))
  list(x = shared_direction, y = mean_strength)
}

# Returns the ordering used
render_for_condition <- function(row_ordering = NULL, column_ordering = NULL) {
  use_covariates <- FALSE
  if(!is.null(row_ordering) & !is.null(column_ordering)) {
    use_covariates <- TRUE
  }
  # Load no-covariate MAP fits
  pattern_str <- "*_bassetfit.rds"
  regexpr_str <- "_bassetfit.rds"
  if(use_covariates) {
    level_dir <- "output/model_fits/ASV_MAP"
  } else {
    level_dir <- "output/model_fits_nocovariates/ASV_MAP"
  }
  model_list <- list.files(path = level_dir, pattern = pattern_str, full.names = TRUE, recursive = FALSE)
  hosts <<- as.vector(sapply(model_list, function(x) { idx <- regexpr(regexpr_str, x); return(substr(x, idx-3, idx-1)) } ))

  Sigmas <- list()
  for(i in 1:length(hosts)) {
    host <- hosts[i]
    cat("Parsing host",i,"/",length(hosts)," (",host,")\n")
    fit <- readRDS(model_list[i])$fit
    fit <- to_clr(fit)
    Sigmas[[host]] <- fit$Sigma[,,1]
  }

  # Build a matrix (hosts x associations) of proportionalities (\rho_p from Quinn et al.)
  # Code pulled from universal_microbes.R

  coordinate_number <- dim(Sigmas[[1]])[1]
  label_pairs <- matrix(NA, coordinate_number, coordinate_number)
  for(i in 1:coordinate_number) {
    for(j in 1:coordinate_number) {
      if(i < j) {
        label_pairs[i,j] <- paste0(i,"_",j)
      }
    }
  }
  label_pairs <- label_pairs[upper.tri(label_pairs, diag=F)]

  associations <- matrix(NA, length(hosts), (coordinate_number^2)/2 - coordinate_number/2)
  for(m in 1:length(hosts)) {
    host <- hosts[m]
    Sigma <- Sigmas[[m]]
    rhos <- c()
    for(j in 2:coordinate_number) {
      for(i in 1:(j-1)) {
        var_i <- Sigma[i,i]
        var_j <- Sigma[j,j]
        var_i_minus_j <- var_i + var_j - 2*Sigma[i,j]
        rho_ij <- 1 - var_i_minus_j / (var_i + var_j)
        rhos <- c(rhos, rho_ij)
      }
    }
    associations[m,] <- rhos
  }

  if(is.null(row_ordering)) {
    # In preparation for clustering, pull the "primary social group" labels we've been using
    data <- load_data(tax_level = "ASV")
    data <- subset_samples(data, sname %in% hosts)
    metadata <- sample_data(data)
    # create a list indexed by host name
    group_labels <- numeric(length(hosts))
    names(group_labels) <- hosts
      
    primary_group <- suppressWarnings(metadata %>%
                                        select(c("sname", "collection_date", "grp")) %>%
                                        filter(sname %in% hosts) %>% 
                                        group_by(sname, grp) %>%
                                        tally() %>%
                                        slice(which.max(n)))
    for(host in hosts) {
      group_labels[host] <- primary_group[primary_group$sname == host,]$grp[[1]]
    }
    group_labels <- as.factor(group_labels)

    # Impose host-ordering by group (clustered within groups)
    groups <- unique(group_labels)
    host.reorder <- c()
    host.reorder.labels <- c()
    for(g in 1:length(groups)) {
      group_idx <- which(group_labels == groups[g])
      d <- dist(associations[group_idx,])
      clustering.hosts <- hclust(d)
      host.reorder.within.group <- group_idx[clustering.hosts$order]
      host.reorder <- c(host.reorder, host.reorder.within.group)
      host.reorder.labels <- c(host.reorder.labels, unname(sapply(hosts[host.reorder.within.group], function(x) {
        paste0(x, " (", group_labels[x], ")")
      })))
    }
    # print(host.reorder.labels)
    associations <- associations[host.reorder,]
  } else {
    associations <- associations[row_ordering,]
  }

  if(is.null(column_ordering)) {
    # Impose association-ordering
    d <- dist(t(associations))
    clustering.assoc <- hclust(d)
    associations <- associations[,clustering.assoc$order]
    label_pairs <- label_pairs[clustering.assoc$order]
  } else {
    associations <- associations[,column_ordering]
    label_pairs <- label_pairs[column_ordering]
  }

  # Render heatmap
  df <- gather_array(associations, "proportionality", "host", "pair")
  p <- ggplot(df, aes(pair, host)) +
    geom_tile(aes(fill = proportionality)) +
    scale_fill_gradient2(low = "darkblue", high = "darkred")
  if(use_covariates) {
    ggsave("rug_cov.png", p, units="in", dpi=150, height=5, width=15)
  } else {
    ggsave("rug_nocov.png", p, units="in", dpi=150, height=5, width=15)
  }

  point_data <- data.frame(x = c(), y = c())
  for(i in 1:ncol(associations)) {
    temp <- calc_xy(associations[,i])
    point_data <- rbind(point_data, data.frame(x = temp$x, y = temp$y))
  }

  bg_data <- data.frame(x = c(), y = c(), z = c())
  shared_direction_proportion <- seq(from = 0.5, to = 1, by = 0.05)
  mean_association_strength <- seq(from = 0, to = 1, by = 0.05)
  for(x in shared_direction_proportion) {
    for(y in mean_association_strength) {
      bg_data <- rbind(bg_data, data.frame(x = x, y = y, z = x*y))
    }
  }

  p <- ggplot() +
    geom_raster(data = bg_data, aes(x = x, y = y, fill = z)) +
    scale_fill_gradient(low = "white", high = "red") +
    geom_point(data = point_data, aes(x = x, y = y), color = "#444444", shape = 1, fill = NA) +
    xlab("proportion shared direction") +
    ylab("mean strength of associations") +
    labs(fill = "weighted \nuniversality\nscore")
  if(use_covariates) {
    ggsave("heatmap_cov.png", p, units = "in", dpi = 150, height = 6, width = 7)
    return(list(row_ordering = host.reorder, column_ordering = clustering.assoc$order))
  } else {
    ggsave("heatmap_nocov.png", p, units = "in", dpi = 150, height = 6, width = 7)
  }
}

orderings <- render_for_condition(row_ordering = NULL, column_ordering = NULL)
render_for_condition(row_ordering = orderings$row_ordering, column_ordering = orderings$column_ordering)

