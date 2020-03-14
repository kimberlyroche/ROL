#' Calculate Riemannian distances between all posterior samples
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param logratio logratio representation to use (e.g. "alr", "ilr", "clr")
#' @return list of MAP covariance estimates indexed by host short name
#' @import driver
#' @import stray
#' @export
#' @examples
#' Sigmas <- load_MAP_estimates(tax_level="genus", logratio="alr")
load_MAP_estimates <- function(tax_level="genus", logratio="alr") {
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

#' Calculate Riemannian distances between all posterior samples
#' 
#' @param taxonomy a named list of taxonomic assignment for a designated sequence variant
#' @param deepest_tax_level the finest resolution taxonomy level to return (if possible)
#' @return readable label of the form "phylum Tenericutes"
#' @export
#' @examples
#' data <- read_file(file.path("input","filtered_genus_5_20.rds"))
#' taxonomy <- tax_table(data)@.Data
#' # get deepest resolve level for the first sequence variant
#' label <- get_deepest_assignment(taxonomy[1,], deepest_tax_level="species")
get_deepest_assignment <- function(taxonomy, deepest_tax_level="genus") {
  for(j in which(names(taxonomy) == deepest_tax_level):1) {
    if(!is.na(taxonomy[j])) {
      return(paste0(names(taxonomy)[j]," ",taxonomy[j]))
    }
  }
  return(NULL)
}

#' Calculate Riemannian distances between all posterior samples
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param logratio logratio representation to use (e.g. "alr", "ilr", "clr")
#' @param other_category "Other" taxonomic category of collapsed low abundance taxa (if available)
#' @details If other_category is provided, these taxa will be rendered as *LR(Other).
#' @return readable list of taxonomic labels (e.g. of the form "CLR(phylum Tenericutes)")
#' @export
#' @examples
#' tax_labels <- assign_concise_taxonomy(tax_level="genus", logratio="alr", other_category=0)
assign_concise_taxonomy <- function(tax_level="genus", logratio="alr", other_category=0) {
  if(logratio != "alr" & logratio != "clr") {
    stop(paste0("Only logratio representations ALR and CLR allowed!\n"))
  }
  data <- read_file(file.path("input",paste0("filtered_",tax_level,"_5_20.rds")))
  taxonomy <- tax_table(data)@.Data
  level_number <- which(colnames(taxonomy) == tax_level)
  if(logratio == "alr") {
    labels <- character(nrow(taxonomy)-1)
    alr.ref <- get_deepest_assignment(taxonomy[nrow(tax),])
  } else {
    labels <- character(nrow(taxonomy))
  }
  for(i in 1:length(labels)) {
    if(i == other_category) {
      if(logratio == "alr") {
        labels[i] <- paste0("ALR(Other / ",alr.ref,")")
      } else {
        labels[i] <- paste0("CLR(Other)")
      }
    } else {
      if(logratio == "alr") {
        labels[i] <- paste0("ALR(",get_deepest_assignment(taxonomy[i,]),")/ALR(",alr.ref,")")
      } else {
        labels[i] <- paste0("CLR(",get_deepest_assignment(taxonomy[i,]),")")
      }
    }
  }
  return(labels)
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
#' interactions <- get_all_vs_one_correlations(taxon_idx=1, tax_level="genus", logratio="alr")
get_all_vs_one_correlations <- function(taxon_idx, tax_level="genus", logratio="alr", Sigmas=NULL) {
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
#' interactions <- get_pairwise_correlations(tax_level="genus", logratio="alr")
get_pairwise_correlations <- function(tax_level="genus", logratio="alr", Sigmas=NULL) {
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

#' Plot heatmap over all pairwise (MAP) correlations between microbes at designated taxonomic level
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param logratio logratio representation to use (e.g. "alr", "ilr", "clr")
#' @param Sigmas optional list (indexed by host short name) of MAP estimates of microbial covariance; if not provided, this will be loaded
#' @param taxon_idx the logratio coordinate to render correlations against; if NULL, render all pairwise correlations
#' @param show_plot show() plot in addition to rendering it to a file
#' @return NULL
#' @import driver
#' @export
#' @examples
#' Sigmas <- load_MAP_estimates(tax_level="genus", logratio="clr")
#' plot_interaction_heatmap(tax_level="genus", logratio="clr", Sigmas=Sigmas)
plot_interaction_heatmap <- function(tax_level="genus", logratio = "alr", Sigmas=NULL, taxon_idx=NULL, show_plot=FALSE) {
  if(logratio != "alr" & logratio != "clr") {
    stop(paste0("Only logratio representations ALR and CLR allowed!\n"))
  }

  if(is.null(taxon_idx)) {
    pairs_obj <- get_pairwise_correlations(tax_level=tax_level, logratio=logratio, Sigmas=Sigmas)
    labels <- pairs_obj$labels
    interactions <- pairs_obj$interactions
    
    # hierarchically cluster all interactions
    d <- dist(t(interactions))
    clustering <- hclust(d)
    # reorder all
    interactions.reordered <- interactions[,clustering$order]
    labels.reordered <- labels[clustering$order]
    df <- gather_array(interactions.reordered, "correlation", "host", "pair")
    # plot
    p <- ggplot(df, aes(pair, host)) +
      geom_tile(aes(fill = correlation)) +
      scale_fill_gradient2(low = "darkblue", high = "darkred")
    if(show_plot) {
      show(p)
    }
    save_dir <- check_output_dir(c("output","plots",paste0(tax_level,"_MAP")))
    ggsave(file.path(save_dir,paste0("microbe_pair_correlations_",logratio,".png")),
           p, units="in", dpi=150, height=5, width=15)
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
#' @param threshold correlation threshold for "interesting" interaction pairs
#' @param show_plot show() plot in addition to rendering it to a file
#' @details Correlations are evaluated in the CLR. A threshold of 0.4 will identify pairs of 
#' taxa with median correlation of > 0.4 or < -0.4. Pairs are rendered as bigraphs.
#' @return NULL
#' @import ggplot2
#' @export
#' @examples
#' get_universal_interactions(tax_level="genus", threshold=0.4)
get_universal_interactions <- function(tax_level="genus", threshold=0.4, show_plot=FALSE) {
  model_list <- get_fitted_model_list(tax_level=tax_level, MAP=TRUE)
  pairs_obj <- get_pairwise_correlations(tax_level=tax_level, logratio="clr")
  labels <- pairs_obj$labels
  interactions <- pairs_obj$interactions
  
  # select compelling interactions by high median absolute correlation
  colMedians <- apply(interactions, 2, function(x) median(x))
  
  criteria <- list(negative=-threshold, positive=threshold)
  for(criterion in names(criteria)) {
    if(criteria[[criterion]] > 0) {
      interesting_pairs <- which(colMedians >= criteria[[criterion]])
    } else {
      interesting_pairs <- which(colMedians <= criteria[[criterion]])
    }
    
    if(length(interesting_pairs) > 0) {
      
      # get readable labels
      tax_labels <- assign_concise_taxonomy(tax_level=tax_level, logratio="clr")
      
      df <- data.frame(x=c(), y=c(), sign=c(), value=c())
      for(p_idx in interesting_pairs) {
        interesting_pair <- labels[p_idx]
        microbe_pair <- as.numeric(strsplit(interesting_pair, "_")[[1]])
        label.1 <- tax_labels[microbe_pair[1]]
        label.2 <- tax_labels[microbe_pair[2]]
        df <- rbind(df, data.frame(x=label.1, y=label.2, sign=sign(colMedians[p_idx]),
                                   value=2*abs(colMedians[p_idx])))
        cat(paste0("Interesting pair (",p_idx,"): ",label.1,", ",label.2,"\n"))
      }
      
      # plot bigraph
      text_df.left <- data.frame(y=1:length(unique(df$x)), label=unique(df$x))
      text_df.right <- data.frame(y=seq(1, length(unique(df$x)), length.out=length(unique(df$y))), label=unique(df$y))
      
      df <- cbind(df, node1=text_df.left$y[as.numeric(df$x)])
      df <- cbind(df, node2=text_df.right$y[as.numeric(df$y)])
      df$sign <- as.factor(df$sign)
      
      if(criterion == "positive") {
        edge_color <- "#379C19"
      } else {
        edge_color <- "#E6194B"
      }
      
      p <- ggplot() +
        geom_segment(data=df, aes(x=1, y=node1, xend=4, yend=node2, color=sign, size=value), alpha=0.2) +
        scale_color_manual(values=c(edge_color)) +
        scale_size_identity() + # use the width specified by `weight`
        geom_text(data=text_df.left, aes(x=1, y=y, label=label), size=4) +
        geom_text(data=text_df.right, aes(x=4, y=y, label=label), size=4) +
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
}
