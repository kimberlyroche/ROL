#' Load MAP estimates of microbial covariances from fitted models at the desired taxonomic level
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param logratio logratio representation to use (e.g. "alr", "ilr", "clr")
#' @return list of MAP covariance estimates indexed by host short name
#' @import driver
#' @import stray
#' @export
#' @examples
#' Sigmas <- load_MAP_estimates(tax_level="ASV", logratio="alr")
load_MAP_estimates <- function(tax_level="ASV", logratio="alr") {
  model_list <- get_fitted_model_list(tax_level=tax_level, MAP=TRUE)
  Sigmas <- list()
  for(i in 1:length(model_list$hosts)) {
    host <- model_list$hosts[i]
    #cat("Processing",host,"\n")
    fit <- readRDS(model_list$model_list[i])$fit
    if(logratio == "clr") {
      fit <- to_clr(fit)
    }
    if(logratio == "ilr") {
      fit <- to_ilr(fit)
    }
    Sigmas[[host]] <- fit$Sigma[,,1]
  }
  return(Sigmas)
}

#' Load full posterior estimates of microbial covariances from fitted models at the desired taxonomic level
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param logratio logratio representation to use (e.g. "alr", "ilr", "clr")
#' @return list of full posterior covariance estimates indexed by host short name
#' @import driver
#' @import stray
#' @export
#' @examples
#' Sigmas <- load_full_posteriors(tax_level="ASV", logratio="alr")
load_full_posteriors <- function(tax_level="ASV", logratio="alr") {
  model_list <- get_fitted_model_list(tax_level=tax_level, MAP=FALSE)
  Sigmas <- list()
  for(i in 1:length(model_list$hosts)) {
    host <- model_list$hosts[i]
    cat("Processing",host,"\n")
    fit <- readRDS(model_list$model_list[i])$fit
    if(logratio == "clr") {
      fit <- to_clr(fit)
    }
    if(logratio == "ilr") {
      fit <- to_ilr(fit)
    }
    Sigmas[[host]] <- fit$Sigma
  }
  return(Sigmas)
}

#' Get correlations between one microbe and all others at designated taxonomic level
#' 
#' @param taxon_idx the logratio coordinate to render correlations against
#' @param tax_level taxonomic level at which to agglomerate data
#' @param logratio logratio representation to use (e.g. "alr", "ilr", "clr")
#' @param Sigmas optional list (indexed by host short name) of MAP estimates of microbial covariance; if not provided, this will be loaded
#' @return correlation matrix of dimensions {number of hosts} x {number of logratio taxa - 1}
#' @export
#' @examples
#' interactions <- get_all_vs_one_correlations(taxon_idx=1, tax_level="ASV", logratio="alr")
get_all_vs_one_correlations <- function(taxon_idx, tax_level="ASV", logratio="alr", Sigmas=NULL) {
  model_list <- get_fitted_model_list(tax_level=tax_level, MAP=TRUE)
  if(logratio == "alr") {
    coordinate_number <- model_list$D-1
  } else {
    coordinate_number <- model_list$D
  }
  if(is.null(Sigmas)) {
    Sigmas <- load_MAP_estimates(tax_level=tax_level, logratio=logratio)
  }
  if(dim(Sigmas[1][[1]])[1] != coordinate_number) {
    stop(paste0("Dimensions of covariance matrices don't match expected logratio dimensions!\n"))
  }
  
  interactions <- matrix(NA, length(model_list$hosts), coordinate_number)
  for(i in 1:length(model_list$hosts)) {
    host <- model_list$hosts[i]
    # convert this host's MAP covariance to correlation
    interactions[i,] <- cov2cor(Sigmas[[host]])[taxon_idx,]
  }
  # omit-self-interactions
  interactions <- interactions[,-taxon_idx]
  return(interactions)
}

#' Get all pairwise (MAP) correlations between microbes at designated taxonomic level
#' 
#' @param taxon_idx the logratio coordinate to render correlations against
#' @param tax_level taxonomic level at which to agglomerate data
#' @param logratio logratio representation to use (e.g. "alr", "ilr", "clr")
#' @param Sigmas optional list (indexed by host short name) of MAP estimates of microbial covariance; if not provided, this will be loaded
#' @return correlation matrix of dimensions {number of hosts} x {number of unique interactions between logratio taxa}
#' @export
#' @examples
#' interactions <- get_pairwise_correlations(tax_level="ASV", logratio="alr")
get_pairwise_correlations <- function(tax_level="ASV", logratio="alr", Sigmas=NULL) {
  model_list <- get_fitted_model_list(tax_level=tax_level, MAP=TRUE)
  if(logratio == "alr") {
    coordinate_number <- model_list$D-1
  } else {
    coordinate_number <- model_list$D
  }
  if(is.null(Sigmas)) {
    Sigmas <- load_MAP_estimates(tax_level=tax_level, logratio=logratio)
  }
  if(dim(Sigmas[1][[1]])[1] != coordinate_number) {
    stop(paste0("Dimensions of covariance matrices don't match expected logratio dimensions!\n"))
  }
  
  # build a map from which we can translate an interaction number with the original pair of taxa
  # this is clumsy but it works
  label_pairs <- matrix(NA, coordinate_number, coordinate_number)
  for(i in 1:coordinate_number) {
    for(j in 1:coordinate_number) {
      if(i < j) {
        label_pairs[i,j] <- paste0(i,"_",j)
      }
    }
  }
  # collect unique interactions only
  label_pairs <- label_pairs[upper.tri(label_pairs, diag=F)]
  
  # plot all in heatmap as hosts x pairs of microbial interaction
  # `.tri` here indicates we're just keeping the upper triangular of the interaction matrix, i.e.
  #   the unique interactions
  interaction_pairs <- matrix(NA, length(model_list$hosts), (coordinate_number^2)/2 - coordinate_number/2)
  for(i in 1:length(model_list$hosts)) {
    host <- model_list$hosts[i]
    # convert this host's MAP covariance to correlation
    corr_mat <- cov2cor(Sigmas[[host]])
    # stack all unique pairwise correlations between microbes in a row associated with this host
    interaction_pairs[i,] <- corr_mat[upper.tri(corr_mat, diag=F)]
  }
  return(list(labels=label_pairs, interactions=interaction_pairs))
}

#' Get all pairwise (MAP) proportionalities between microbes at designated taxonomic level
#' 
#' @param taxon_idx the logratio coordinate to render correlations against
#' @param tax_level taxonomic level at which to agglomerate data
#' @param logratio logratio representation to use (e.g. "alr", "ilr", "clr")
#' @param Sigmas optional list (indexed by host short name) of MAP estimates of microbial covariance; if not provided, this will be loaded
#' @return correlation matrix of dimensions {number of hosts} x {number of unique interactions between logratio taxa}
#' @export
#' @examples
#' interactions <- get_pairwise_proportionalities(tax_level = "ASV", logratio = "alr")
get_pairwise_proportionalities <- function(tax_level = "ASV", logratio = "alr", Sigmas = NULL) {
  model_list <- get_fitted_model_list(tax_level = tax_level, MAP = TRUE)
  if(logratio == "alr") {
    coordinate_number <- model_list$D-1
  } else {
    coordinate_number <- model_list$D
  }
  if(is.null(Sigmas)) {
    Sigmas <- load_MAP_estimates(tax_level=tax_level, logratio=logratio)
  }
  if(dim(Sigmas[1][[1]])[1] != coordinate_number) {
    stop(paste0("Dimensions of covariance matrices don't match expected logratio dimensions!\n"))
  }
  
  # build a map from which we can translate an interaction number with the original pair of taxa
  # this is clumsy but it works
  label_pairs <- matrix(NA, coordinate_number, coordinate_number)
  for(i in 1:coordinate_number) {
    for(j in 1:coordinate_number) {
      if(i < j) {
        label_pairs[i,j] <- paste0(i,"_",j)
      }
    }
  }
  # collect unique interactions only
  label_pairs <- label_pairs[upper.tri(label_pairs, diag=F)]
  
  # plot all in heatmap as hosts x pairs of microbial interaction
  # `.tri` here indicates we're just keeping the upper triangular of the interaction matrix, i.e.
  #   the unique interactions
  interaction_pairs <- matrix(NA, length(model_list$hosts), (coordinate_number^2)/2 - coordinate_number/2)
  for(m in 1:length(model_list$hosts)) {
    host <- model_list$hosts[m]
    # convert this host's MAP covariance to correlation
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
    # stack all unique pairwise correlations between microbes in a row associated with this host
    interaction_pairs[m,] <- rhos
  }
  return(list(labels=label_pairs, interactions=interaction_pairs))
}

#' Plot heatmap over all pairwise (MAP) correlations between microbes at designated taxonomic level
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param logratio logratio representation to use (e.g. "alr", "ilr", "clr")
#' @param Sigmas optional list (indexed by host short name) of MAP estimates of microbial covariance; if not provided, this will be loaded
#' @param cluster optional flag to hierarchically cluster across hosts and interactions
#' @param taxon_idx optional logratio coordinate to render correlations against; if NULL, render all pairwise correlations
#' @param show_plot optional flag to show() plot in addition to rendering it to a file
#' @param return_matrix optional flag to return host x interaction correlation value matrix
#' @return NULL or heatmap matrix
#' @import driver
#' @export
#' @examples
#' Sigmas <- load_MAP_estimates(tax_level="ASV", logratio="clr")
#' plot_interaction_heatmap(tax_level="ASV", logratio="clr", Sigmas=Sigmas)
plot_interaction_heatmap <- function(tax_level="ASV", logratio = "alr", Sigmas=NULL,
                                     cluster=TRUE, taxon_idx=NULL, show_plot=FALSE, return_matrix=FALSE) {
  if(logratio != "alr" & logratio != "clr") {
    stop(paste0("Only logratio representations ALR and CLR allowed!\n"))
  }

  if(is.null(taxon_idx)) {
    # CLR correlation
    # pairs_obj <- get_pairwise_correlations(tax_level = tax_level, logratio = logratio, Sigmas = Sigmas)
    # proportionality (\rho_p as in Quinn et a.)
    pairs_obj.prop <- get_pairwise_proportionalities(tax_level = tax_level, logratio = logratio, Sigmas = Sigmas)
    labels <- pairs_obj$labels
    interactions <- pairs_obj$interactions
    
    # hierarchically cluster all interactions
    if(cluster) {
      d <- dist(interactions)
      clustering.hosts <- hclust(d)
      d <- dist(t(interactions))
      clustering.interactions <- hclust(d)
      # reorder all
      interactions.reordered <- interactions[clustering.hosts$order,]
      interactions.reordered <- interactions.reordered[,clustering.interactions$order]
      labels.reordered <- labels[clustering.interactions$order]
    } else {
      interactions.reordered <- interactions
      labels.reordered <- labels
    }
    df <- gather_array(interactions.reordered, "proportionality", "host", "pair")
    # plot
    p <- ggplot(df, aes(pair, host)) +
      geom_tile(aes(fill = proportionality)) +
      scale_fill_gradient2(low = "darkblue", high = "darkred")
    if(show_plot) {
      show(p)
    }
    save_dir <- check_output_dir(c("output","plots",paste0(tax_level,"_MAP")))
    ggsave(file.path(save_dir,paste0("microbe_pair_proportionalities_",logratio,".png")),
           p, units="in", dpi=150, height=5, width=15)
        
    if(return_matrix) {
      return(list(clustering.hosts = clustering.hosts,
                  clustering.interactions = clustering.interactions,
                  interactions.reordered = interactions.reordered,
                  labels.reordered = labels.reordered))
    }
  } else {
    interactions <- get_all_vs_one_correlations(taxon_idx, tax_level=tax_level, logratio=logratio, Sigmas=Sigmas)
    
    # plot interactions in order
    df <- gather_array(interactions, "correlation", "host", "pair")
    p <- ggplot(df, aes(pair, host)) +
      geom_tile(aes(fill = correlation)) +
      scale_fill_gradient2(low = "darkblue", high = "darkred")
    if(show_plot) {
      show(p)
    }
    save_dir <- check_output_dir(c("output","plots",paste0(tax_level,"_MAP")))
    ggsave(file.path(save_dir,paste0("microbe_",taxon_idx,"_correlations_",logratio,".png")),
           p, units="in", dpi=150, height=5, width=6)
  }
}

#' Identify and plot "universal" interactions at the designated taxonomic level
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param show_plot show() plot in addition to rendering it to a file
#' @param order_by ordering of bigraph components; options are "abundance" or "taxonomy"
#' @details Correlations are evaluated in the CLR. A threshold of 0.4 will identify pairs of 
#' taxa with median correlation of > 0.4 or < -0.4. Pairs are rendered as bigraphs.
#' @return data.frame of microbial pairs (by index in taxonomy and label) and median correlation across hosts
#' @import ggplot2
#' @import phyloseq
#' @export
#' @examples
#' interaction_matrix <- get_universal_interactions(tax_level="ASV")
get_universal_interactions <- function(tax_level="ASV", show_plot=FALSE, order_by="abundance") {
  if(order_by != "abundance" & order_by != "taxonomy") {
    stop(paste0("Invalid interaction ordering '",order_by,"'!\n"))
  }
  model_list <- get_fitted_model_list(tax_level=tax_level, MAP=TRUE)
  pairs_obj <- get_pairwise_correlations(tax_level=tax_level, logratio="clr")
  labels <- pairs_obj$labels
  interactions <- pairs_obj$interactions
  
  # select compelling interactions by high median absolute correlation
  colMedians <- apply(interactions, 2, function(x) median(x))

  # order the interactions by median (large - to large +); we'll pull out the strongest
  ranks <- order(colMedians)

  # get average abundance; we'll order interactions (x vs. y) by this
  data <- load_data(tax_level=tax_level)
  if(order_by == "abundance") {
    avg_abundance <- colMeans(otu_table(data)@.Data)
    names(avg_abundance) <- NULL
  } else {
    alr_ref <- formalize_parameters(data)$alr_ref
    taxonomy <- get_taxonomy(data, alr_ref)
    colnames(taxonomy) <- NULL
    taxonomy[is.na(taxonomy)] <- "zzz"
    tax_names <- apply(taxonomy, 1, function(x) paste(x, collapse="/"))
  }

  criteria <- c("negative", "positive")
  for(criterion in criteria) {
    # take the top and bottom 10% of interactions
    top_k <- min(round(0.1*length(ranks)), 100)
    if(criterion == "positive") {
      # positive interactions
      interesting_pairs <- ranks[(length(ranks)-top_k+1):length(ranks)]
      # cat(paste0("Evaluating top ",top_k," positive interactions...\n"))
    } else {
      # negative interactions
      interesting_pairs <- ranks[1:top_k]
      # cat(paste0("Evaluating top ",top_k," negative interactions...\n"))
    }
    
    if(length(interesting_pairs) > 0) {
      # get readable labels
      tax_labels <- assign_concise_taxonomy(tax_level=tax_level, logratio="clr")
      
      df <- data.frame(x=c(), x_idx=c(), y=c(), y_idx=c(), sign=c(), value=c())
      for(p_idx in interesting_pairs) {
        interesting_pair <- labels[p_idx]
        microbe_pair <- as.numeric(strsplit(interesting_pair, "_")[[1]])
        if(tax_level == "ASV") {
          label.1 <- paste0("Microbe #",microbe_pair[1]," in ",tax_labels[microbe_pair[1]])
          label.2 <- paste0("Microbe #",microbe_pair[2]," in ",tax_labels[microbe_pair[2]])
        } else {
          label.1 <- tax_labels[microbe_pair[1]]
          label.2 <- tax_labels[microbe_pair[2]]
        }
        df <- rbind(df, data.frame(x=label.1, x_idx=microbe_pair[1], y=label.2, y_idx=microbe_pair[2],
                                   sign=sign(colMedians[p_idx]), value=2*abs(colMedians[p_idx])))
        # cat(paste0("Interesting pair (",p_idx,"): ",label.1,", ",label.2,"\n"))
      }

      x_nodes <- unique(df$x_idx)
      x_labels <- unique(df$x)
      y_nodes <- unique(df$y_idx)
      y_labels <- unique(df$y)
      if(order_by == "abundance") {
        # get order of integer labels and text labels according to average abundance
        # for x
        x_ab <- avg_abundance[x_nodes]
        ab_order <- order(x_ab, decreasing=FALSE)
        x_ordered <- x_nodes[ab_order]
        x_labels_ordered <- x_labels[ab_order]
        # for y
        y_ab <- avg_abundance[y_nodes]
        ab_order <- order(y_ab, decreasing=FALSE)
        y_ordered <- y_nodes[ab_order]
        y_labels_ordered <- y_labels[ab_order]
      } else {
        # get order of integer labels and text labels according to taxonomy
        # for x
        x_tax <- tax_names[x_nodes]
        tax_order <- order(x_tax, decreasing=TRUE)
        x_ordered <- x_nodes[tax_order]
        x_labels_ordered <- x_labels[tax_order]
        # for y
        y_tax <- tax_names[y_nodes]
        tax_order <- order(y_tax, decreasing=TRUE)
        y_ordered <- y_nodes[tax_order]
        y_labels_ordered <- y_labels[tax_order]
      }

      # apply the ordering
      text_df.left <- data.frame(ypos_left=1:length(x_ordered),
                                 x_idx=x_ordered,
                                 label_left=x_labels_ordered)
      text_df.right <- data.frame(ypos_right=seq(1, length(x_ordered), length.out=length(y_ordered)),
                                  y_idx=y_ordered,
                                  label_right=y_labels_ordered)
      
      # join the interactions in 'df' and the nodes in 'text_df.*'
      df <- merge(df, text_df.left, by="x_idx")
      df <- merge(df, text_df.right, by="y_idx")
      df$sign <- as.factor(df$sign)
      
      if(criterion == "positive") {
        edge_color <- "#E6194B"
      } else {
        edge_color <- "#202ABA"
      }
      
      p <- ggplot() +
        geom_segment(data=df, aes(x=1, y=ypos_left, xend=4, yend=ypos_right, color=sign, size=value), alpha=0.2) +
        scale_color_manual(values=c(edge_color)) +
        scale_size_identity() + # use the width specified by `weight`
        geom_text(data=text_df.left, aes(x=1, y=ypos_left, label=label_left), size=4) +
        geom_text(data=text_df.right, aes(x=4, y=ypos_right, label=label_right), size=4) +
        xlim(0, 5) +
        theme(legend.position = "none",
              panel.grid = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              panel.background = element_blank()) 
      if(show_plot) {
        show(p)
      }
      save_dir <- check_output_dir(c("output","plots",paste0(tax_level,"_MAP")))
      ggsave(file.path(save_dir,paste0(criterion,"_interactions_",tax_level,".png")),
             p, units="in", dpi=150, height=10, width=10)

    }
  }

  # messy but assemble the data.frame to return
  interaction_set <- as.data.frame(t(sapply(pairs_obj$labels, function(x) as.numeric(strsplit(x, "_")[[1]]))))
  label_obj <- as.data.frame(sapply(interaction_set[,1], function(x) tax_labels[x]))
  label_obj <- cbind(label_obj, sapply(interaction_set[,2], function(x) tax_labels[x]))
  interaction_set <- cbind(interaction_set, label_obj, colMedians)
  rownames(interaction_set) <- NULL
  colnames(interaction_set) <- c("idx.1", "idx.2", "label.1", "label.2", "median_corr")
  # reorder by ranks
  interaction_set <- interaction_set[ranks,]
  return(interaction_set)
}

# old version
calc_universality_score.old <- function(x) {
  max_sd <- sd(c(rep(-1, length(x)/2), rep(1, length(x)/2)))
  return(((max_sd - sd(x))/max_sd) * abs(mean(x)))
}

#' Calculate a score (0.0 to 1.0) for the "universality" of given microbial interaction
#' 
#' @param x vector of correlations between a pair of CLR microbes across individuals
#' @details An element of x should be the (MAP) correlation between a pair of CLR microbes.
#' The length of x is assumed equal to the number of hosts.
#' @return score
#' @export
#' @examples
#' interaction_matrix <- get_universal_interactions(tax_level="ASV")
#' score <- calc_universality_score(interaction_matrix[,1])
calc_universality_score <- function(x) {
  x.sign <- sapply(x, sign)
  neg.idx <- which(x.sign < 0)
  pos.idx <- which(x.sign > 0)
  neg.abs.median <- NA
  pos.abs.median <- NA
  if(length(neg.idx) > 0) {
    neg.abs.median <- abs(median(x[neg.idx]))
  }
  if(length(pos.idx) > 0) {
    pos.abs.median <- abs(median(x[pos.idx]))
  }
  score <- 0
  # R doesn't do short circuit evaluation? yes, it does: update with &&
  if(!is.na(neg.abs.median) & !is.na(pos.abs.median)) {
    if(pos.abs.median > neg.abs.median) {
      # use this as majority direction
      score <- (length(pos.idx)/length(x))*(pos.abs.median)
    } else {
      score <- (length(neg.idx)/length(x))*(neg.abs.median)
    }
  } else {
    if(is.na(neg.abs.median)) {
      score <- (length(pos.idx)/length(x))*(pos.abs.median)
    } else {
      score <- (length(neg.idx)/length(x))*(neg.abs.median)
    }
  }
  score
}

#' Calculate a score (0.0 to 1.0) for the "universality" of given microbe
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param Sigmas optional list (indexed by host short name) of MAP estimates of microbial covariance; if not provided, this will be loaded
#' @details These scores do not currently incorporate posterior uncertainty. Returns results for correlations in CLR.
#' @return named list of per-microbe scores
#' @export
#' @examples
#' Sigmas <- load_MAP_estimates(tax_level="ASV", logratio="clr")
#' scores <- calc_microbe_universality_score(tax_level="ASV", Sigmas=Sigmas)
#' reorder <- order(unlist(scores), decreasing=TRUE)
#' scores <- unlist(scores)[reorder]
calc_microbe_universality_score <- function(tax_level="ASV", Sigmas=NULL) {
  if(is.null(Sigmas)) {
    Sigmas <- load_MAP_estimates(tax_level=tax_level, logratio="clr")
  }
  n_taxa <- nrow(Sigmas[[1]])
  hosts <- names(Sigmas)
  # transform these posterior covariance estimates to correlation
  Sigmas_corr <- Sigmas
  for(host in hosts) {
    Sigmas_corr[[host]] <- cov2cor(Sigmas_corr[[host]])
  }

  # get labels
  data <- load_data(tax_level=tax_level)
  alr_ref <- formalize_parameters(data)$alr_ref
  tax <- get_taxonomy(data, alr_ref)

  # get scores
  median_scores <- list()
  for(focal_microbe in 1:n_taxa) {
    focal_scores <- c()
    for(other_microbe in setdiff(1:n_taxa, focal_microbe)) {
      focal_interactions <- c()
      for(host in hosts) {
        focal_interactions <- c(focal_interactions, Sigmas_corr[[host]][focal_microbe, other_microbe])
      }
      focal_scores <- c(focal_scores, calc_universality_score(focal_interactions))
    }
    label <- get_deepest_assignment(tax[focal_microbe,], deepest_tax_level="genus")
    median_scores[[label]] <- median(focal_scores)
  }
  return(median_scores)
}

