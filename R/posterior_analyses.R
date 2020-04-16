#' Calculate distances between all posterior samples of either Sigma or Lambda
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param which_measure estimated object to embed, either Lambda or Sigma
#' @param which_distance ignored for Lambda; for Sigma allowable options are "Riemannian" (default) or "Euclidean"
#' @param MAP use MAP estimate model output instead of full posterior output
#' @param spike_in copy and permute the samples and calculate distances to these too, so as to give upper bound on distances
#' @details Distance matrix between posterior samples saved to designated output directory. Saves
#' a list containing distance matrix and labels of each row or column (by host).
#' @return NULL
#' @import driver
#' @import stray
#' @import matrixsampling
#' @export
#' @examples
#' calc_posterior_distances(tax_level="ASV", which_measure="Sigma", MAP=FALSE)
calc_posterior_distances <- function(tax_level="ASV", which_measure="Sigma",
  which_distance="Riemannian", MAP=FALSE, spike_in=FALSE) {
  # grab all fitted models
  model_list <- get_fitted_model_list(tax_level=tax_level, MAP=MAP)
  P <- model_list$D - 1 # ALR
  n_samples <- model_list$n_samples
  n_hosts <- length(model_list$hosts)
  # initialize samples matrix
  if(which_measure == "Sigma") {
    if(which_distance == "Riemannian") {
      if(MAP & spike_in) {
        # make room for 2x the samples
        all_samples <- matrix(NA, P, P*n_samples*2*n_hosts)
      } else {
        all_samples <- matrix(NA, P, P*n_samples*n_hosts)
      }
    } else {
      all_samples <- matrix(NA, (P^2)/2 + P/2, n_samples*n_hosts)
    }
  } else {
    # we'll use per-sample average to manage individuals having different N so each individual's
    # posterior will be summarized as one mean vector
    all_samples <- matrix(NA, n_hosts*n_samples, P)
  }
  # insert samples (column-wise) into samples matrix
  host_labels <- c()
  for(i in 1:n_hosts) {
    fit <- read_file(model_list$model_list[i])$fit
    # convert to ILR; this can be removed
    # V <- driver::create_default_ilr_base(ncategories(fit))
    # fit.ilr <- to_ilr(fit, V)
    # Lambda <- fit.ilr$Lambda
    # Sigma <- fit.ilr$Sigma
    Lambda <- fit$Lambda
    Sigma <- fit$Sigma
    if(which_measure == "Sigma") {
      if(which_distance == "Riemannian") {
        Sigma <- Sigma[,,1:n_samples]
        # symmetrize this guy
        Sigma <- (Sigma + t(Sigma))/2
        if(MAP & spike_in) {
          # here we know n_samples == 1
          offset1 <- (i-1)*P*2+1
          offset2 <- offset1 + P - 1
          all_samples[,offset1:offset2] <- Sigma
          upsilon <- P + 2
          random_Sigma <- rinvwishart(1, upsilon, Sigma*(upsilon - P - 1))[,,1]
          offset3 <- offset2 + 1
          offset4 <- offset3 + P - 1
          all_samples[,offset3:offset4] <- random_Sigma
        } else {
          all_samples[,((i-1)*P*n_samples+1):(i*P*n_samples)] <- Sigma
        }
      } else {
        host_offset <- (i-1)*n_samples
        for(j in 1:n_samples) {
          Sigma_sample <- Sigma[,,j]
          all_samples[,(host_offset+j)] <- c(Sigma_sample[upper.tri(Sigma_sample, diag=T)])
        }
      }
      if(MAP & spike_in) {
        host_labels <- c(host_labels, c(model_list$hosts[i], "random"))
      } else {
        host_labels <- c(host_labels, rep(model_list$hosts[i], n_samples))
      }
    } else {
      collLambda <- t(apply(Lambda, 3, function(X) { apply(X, 1, mean) })) # n_samples x P
      all_samples[((i-1)*n_samples+1):(i*n_samples),] <- collLambda
      host_labels <- c(host_labels, rep(model_list$hosts[i], n_samples))
    }
  }
  
  save_dir <- check_output_dir(c("output"))
  if(MAP) {
    if(spike_in) {
      dist_filename <- file.path(save_dir,paste0(which_measure,"_distance_",tax_level,"_MAP_spikein.rds"))
    } else {
      dist_filename <- file.path(save_dir,paste0(which_measure,"_distance_",tax_level,"_MAP.rds"))
    }
  } else {
    dist_filename <- file.path(save_dir,paste0(which_measure,"_distance_",tax_level,".rds"))
  }
  if(which_measure == "Sigma") {
    if(which_distance == "Riemannian") {
      if(MAP & spike_in) {
        distance_mat <- Riemann_dist_samples(all_samples, n_hosts, n_samples*2)
      } else {
        distance_mat <- Riemann_dist_samples(all_samples, n_hosts, n_samples)
      }
    } else {
      distance_mat <- as.matrix(dist(t(all_samples)))
    }
  } else {
    # use Euclidean distance
    distance_mat <- as.matrix(dist(all_samples))
  }
  saveRDS(list(host_labels=host_labels, distance_mat=distance_mat), file=dist_filename)
}

#' Calculate Riemannian distances between one sample and all others; this is to allow extremely
#' large distance matrices to be calculated in blocks by independent parallel jobs on the
#' cluster
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param sample_idx the posterior sample to be considered against all others
#' @details Output is saved to a file.
#' @return NULL
#' @import driver
#' @import stray
#' @export
#' @examples
#' calc_posterior_distances_row(tax_level="ASV", sample_idx=1)
calc_posterior_distance_row <- function(tax_level="ASV", sample_idx=1) {
  # grab all fitted models
  model_list <- get_fitted_model_list(tax_level=tax_level, MAP=FALSE)
  P <- model_list$D - 1 # ALR
  n_samples <- model_list$n_samples
  if(sample_idx > n_samples) {
    sample_idx <- 1
  }
  n_hosts <- length(model_list$hosts)
  # initialize samples matrix
  # for reference samples, we just want to consider sample indices > this index since
  # we can copy the upper triangular part of the distance matrix into the bottom triangle
  # and save work
  reference_samples <- matrix(NA, P, P*n_hosts)
  all_samples <- matrix(NA, P, P*n_samples*n_hosts)
  # insert samples (column-wise) into samples matrix
  for(i in 1:n_hosts) {
    fit <- read_file(model_list$model_list[i])$fit
    # convert to ILR; this can be removed
    V <- driver::create_default_ilr_base(ncategories(fit))
    fit.ilr <- to_ilr(fit, V)
    Sigma <- fit.ilr$Sigma
    Sigma_ref <- Sigma[,,sample_idx]
    Sigma_full <- Sigma[,,1:n_samples]
    offset1 <- (i-1)*P + 1
    offset2 <- offset1 + P - 1
    reference_samples[,offset1:offset2] <- Sigma_ref
    for(j in 1:n_samples) {
      offset1 <- (j-1)*n_hosts*P + ((i-1)*P) + 1
      offset2 <- offset1 + P - 1
      all_samples[,offset1:offset2] <- Sigma_full[,,j]
    }
  }
  save_dir <- check_output_dir(c("output"))
  dist_filename <- file.path(save_dir,paste0("Sigma_distance_",tax_level,"_",sample_idx,".rds"))
  distance_mat <- Riemann_dist_sets(reference_samples, all_samples, n_hosts, 1, n_samples)
  saveRDS(list(host_labels=model_list$hosts, distance_mat=distance_mat), file=dist_filename)
}

#' Embeds posterior samples using MDS and a pre-calculated distance matrix
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param which_measure estimated object to embed, either Lambda or Sigma
#' @param MAP use MAP estimate model output instead of full posterior output
#' @param spike_in copy and permute the samples and calculate distances to these too, so as to give upper bound on distances
#' @details Distance matrix between posterior samples must be present in designated output directory
#' @return NULL
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @export
#' @examples
#' embed_posteriors(tax_level="ASV", which_measure="Sigma", MAP=FALSE)
embed_posteriors <- function(tax_level="ASV", which_measure="Sigma", MAP=FALSE, spike_in=FALSE) {
  if(MAP) {
    if(spike_in) {
      dist_filename <- file.path("output",paste0(which_measure,"_distance_",tax_level,"_MAP_spikein.rds"))
    } else {
      dist_filename <- file.path("output",paste0(which_measure,"_distance_",tax_level,"_MAP.rds"))
    }
  } else {
    dist_filename <- file.path("output",paste0(which_measure,"_distance_",tax_level,".rds"))
  }
  if(!file.exists(dist_filename)) {
    stop(paste0("Error: unable to locate distance matrix over posterior samples at level ",tax_level,"!\n"))
  }
  dist_obj <- read_file(dist_filename)
  host_labels <- dist_obj$host_labels
  distance_mat <- dist_obj$distance_mat
  
  k <- 6
  if(nrow(distance_mat) < k) {
    k <- nrow(distance_mat)
  }
  fit <- cmdscale(distance_mat, eig=TRUE, k=k-1)
  # I believe in this case the magnitude of the eigenvalues is proportional to the variance
  # explained (technically it differs by a factor of n-1 [the DOF] I think)
  eig_tot <- sum(abs(fit$eig))
  for(i in 1:(k-1)) {
    cat(paste0("Eigenvalue #",i,": ",round(fit$eig[i],2)," (% variance: ",round(abs(fit$eig[i])/eig_tot,2),")\n"))
  }

  # save coordinates of interest in a data.frame; this will have the form
  #   coord     value labels
  # 1     1 -1.876277    ZIB
  # 2     1 -2.336495    ZIB
  # ...
  df <- data.frame(coord=c(), value=c(), labels=c())
  for(i in 1:ncol(fit$points)) {
    df <- rbind(df, data.frame(coord=rep(i, nrow(fit$points)), value=fit$points[,i], labels=host_labels))
  }
  if(MAP) {
    save_dir <- check_output_dir(c("output","plots",paste0(tax_level,"_MAP")))
    if(spike_in) {
      saveRDS(df, file.path("output","plots",paste0(tax_level,"_MAP"),paste0(which_measure,"_ordination_spikein.rds")))
    } else {
      saveRDS(df, file.path("output","plots",paste0(tax_level,"_MAP"),paste0(which_measure,"_ordination.rds")))
    }
  } else {
    save_dir <- check_output_dir(c("output","plots",tax_level))
    saveRDS(df, file.path("output","plots",tax_level,paste0(which_measure,"_ordination.rds")))
  }

  df_centroids <- NULL
  for(i in 1:max(df$coord)) {
    temp <- df[df$coord == i,] %>%
      group_by(labels) %>%
      summarise(mean=mean(value))
    names(temp) <- c("labels", i)
    if(is.null(df_centroids)) {
      df_centroids <- temp
    } else {
      df_centroids <- left_join(df_centroids, temp, by="labels")
    }
  }
  df_centroids <- as.data.frame(df_centroids)
  # transform these such that they're in the same format as the coordinates data.frame
  df_centroids <- gather(df_centroids, "coord", "value", 2:ncol(df_centroids))
  df_centroids <- df_centroids[,colnames(df)]
  
  if(MAP) {
    if(spike_in) {
      saveRDS(df_centroids, file.path("output","plots",paste0(tax_level,"_MAP"),paste0(which_measure,"_ordination_centroids_spikein.rds")))
    } else {
      saveRDS(df_centroids, file.path("output","plots",paste0(tax_level,"_MAP"),paste0(which_measure,"_ordination_centroids.rds")))
    }
  } else {
    saveRDS(df_centroids, file.path("output","plots",tax_level,paste0(which_measure,"_ordination_centroids.rds")))
  }
}

#' Loads fitness annotations associated with hosts from file 
#' 
#' @param hosts a vector of hosts by sname
#' @details Output is a data.frame of fitness annotations sorted by host short name
#' @return data.frame
#' @export
#' @examples
#' outcomes <- load_outcomes(data)
load_outcomes <- function(hosts) {
  outcomes <- read.csv(file.path("input","individual_traits.csv"), header=TRUE)
  outcomes <- outcomes[outcomes$sname %in% hosts,]
  # filter to NA-less measures
  outcomes <- outcomes[,apply(outcomes, 2, function(x) sum(is.na(x))==0)]
  return(outcomes)
}

#' Get other annotations (e.g. sample number) for hosts in data set
#' 
#' @param centroids data.frame of per-host centroids from posterior embedding
#' @param tax_level taxonomic level at which to agglomerate data
#' @param annotation label to assign (e.g. individual)
#' @param MAP use MAP estimate model output instead of full posterior output
#' @details Output is the centroids data.frame appended with the desired annotation. Available annotations are:
#' group, matgroup, samplenumber, sampledensity, mom, dad, momrank, drought, largegroup, 
#' momdied, competingsib, earlyadversity, birthrate_all, birthrate_surviving, alphadiv,
#' wetproportion 
#' @return data.frame
#' @import phyloseq
#' @import dplyr
#' @export
#' @examples
#' centroids <- readRDS(file.path("output","plots","genus","Sigma_ordination_centroids.rds"))
#' labelled_centoids <- get_other_labels(centroids, tax_level="ASV", annotation="group", MAP=FALSE)
get_other_labels <- function(centroids, tax_level="ASV", annotation="group", MAP=FALSE) {
  # pull unique hosts from centroids
  hosts <<- unique(as.character(centroids$labels))
  data <- load_data(tax_level)
  data <- subset_samples(data, sname %in% hosts)
  metadata <- sample_data(data)
  # create a list indexed by host name
  labels <- numeric(length(hosts))
  names(labels) <- hosts
  
  if(annotation == "group") {
    primary_group <- suppressWarnings(metadata %>%
                                        select(c("sname", "collection_date", "grp")) %>%
                                        filter(sname %in% hosts) %>% 
                                        group_by(sname, grp) %>%
                                        tally() %>%
                                        slice(which.max(n)))
    for(host in hosts) {
      labels[host] <- primary_group[primary_group$sname == host,]$grp[[1]]
    }
    labels <- as.factor(labels)
  }
  if(annotation == "matgroup") {
    for(host in hosts) {
      labels[centroids$labels == host] <- metadata[metadata$sname == host,]$matgrp[[1]]
    }
    labels <- as.factor(labels)
  }
  if(annotation == "samplenumber" | annotation == "sampledensity") {
    for(host in hosts) {
      host <<- host
      host_subset <- subset_samples(data, sname == host)
      sample_count <- phyloseq::nsamples(host_subset)
      labels[centroids$labels == host] <- round(sample_count, -1) # discretize
    }
  }
  if(annotation == "sampledensity") {
    md_subset <- suppressWarnings(metadata %>%
                                    select(c("sname", "collection_date")) %>%
                                    filter(sname %in% hosts) %>%
                                    group_by(sname) %>%
                                    summarize(delta=difftime(max(collection_date), min(collection_date), units="days")))
    for(host in hosts) {
      labels[names(labels) == host] <- labels[names(labels) == host]/md_subset[md_subset$sname == host,]$delta[[1]]
    }
    labels <- round(labels*100)
  }
  if(annotation %in% c("mom", "dad")) {
    outcomes <- load_outcomes(hosts)
    for(host in hosts) {
      if(annotation == "mom") { label <- as.character(outcomes[outcomes$sname == host,]$mom) }
      else { label <- as.character(outcomes[outcomes$sname == host,]$dad) }
      if(label == "") { label <- NA }
      labels[names(labels) == host] <- label
    }
  }
  if(annotation %in% c("momrank", "drought", "largegroup", "momdied", "competingsib", "earlyadversity")) {
    outcomes <- load_outcomes(hosts)
    for(host in hosts) {
      if(annotation == "momrank") { labels[names(labels) == host] <- outcomes[outcomes$sname == host,]$mom_lowQuartRank }
      if(annotation == "drought") { labels[names(labels) == host] <- outcomes[outcomes$sname == host,]$bornInDrought }
      if(annotation == "largegroup") { labels[names(labels) == host] <- outcomes[outcomes$sname == host,]$born_largeGroup }
      if(annotation == "momdied") { labels[names(labels) == host] <- outcomes[outcomes$sname == host,]$mom_died }
      if(annotation == "competingsib") { labels[names(labels) == host] <- outcomes[outcomes$sname == host,]$has_CompetingSib }
      if(annotation == "earlyadversity") { labels[names(labels) == host] <- outcomes[outcomes$sname == host,]$EarlyAdversityScore }
    }
    labels <- as.factor(labels)
  }
  if(annotation %in% c("birthrate_all", "birthrate_surviving")) {
    outcomes <- load_outcomes(hosts)
    years <- sapply(as.vector(outcomes$birth_date), function(x) {
      year <- as.numeric(strsplit(x, "/")[[1]][3]);
      if(year < 10) { year <- year + 2000 }
      else { year <- year + 1900 }
    })
    names(years) <- outcomes$sname
    
    for(host in hosts) {
      if(outcomes[outcomes$sname==host,]$sex == "F") {
        years_obs_cont <- difftime(max(metadata[metadata$sname==host,]$collection_date), min(metadata[metadata$sname==host,]$collection_date), units="weeks")/52
        if(years_obs_cont >= 1) {
          if(annotation == "birthrate_all") {
            births <- outcomes[outcomes$sname==host,]$num_live_births_RAW
          } else {
            births <- outcomes[outcomes$sname==host,]$num_surv_births_RAW
          }
          score <- births/as.numeric(years_obs_cont)
          # labels[names(labels) == host] <- round(score,1)
          labels[names(labels) == host] <- round(score*2)/2 # round to nearest half
        } else {
          labels[names(labels) == host] <- NA
        }
      } else {
        labels[names(labels) == host] <- NA
      }
    }
    labels <- as.factor(labels)
  }
  if(annotation == "alphadiv") {
    for(host in hosts) {
      host <<- host
      host_data <- subset_samples(data, sname == host)
      counts <- otu_table(host_data)@.Data # samples x taxa
      props <- apply(counts, 1, function(x) x/sum(x)) # flipped after this: taxa x samples
      props <- props + min(props[props != 0])*0.1 # pseudocount 10x smaller than non-zero min
      alpha_div_shannon <- apply(props, 2, function(x) -sum(x*log2(x)))
      labels[centroids$labels == host] <- round(mean(alpha_div_shannon), 1) # discretize
    }
    labels <- as.factor(labels)
  }
  if(annotation == "wetproportion") {
    for(host in hosts) {
      host <<- host
      host_data <- subset_samples(data, sname == host)
      season_vec <- metadata$season
      labels[centroids$labels == host] <- round(sum(season_vec == "Wet")/length(season_vec),1) # discretize
    }
  }
  
  centroids$labels <- labels
  return(centroids)
}

#' Plot a pair of principal coordinates
#' 
#' @param coordinates data.frame of coordinates from posterior embedding
#' @param centroids data.frame of centroids from posterior embedding
#' @param tax_level taxonomic level at which to agglomerate data
#' @param axis1 PCoA coordinate to display on x-axis
#' @param axis2 PCoA coordinate to display on y-axis
#' @param annotation label to assign (e.g. individual)
#' @param MAP use MAP estimate model output instead of full posterior output
#' @param show_plot show() plot in addition to rendering it to a file
#' @details Coordinates parameter is only necessary where we're plotting full host posteriors (i.e. where
#' annotation="host").
#' @return NULL
#' @import ggplot2
#' @export
#' @examples
#' tax_level <- "genus"
#' annotation <- "group"
#' MAP <- FALSE
#' centroids <- readRDS(file.path("output","plots",tax_level,"Sigma_ordination_centroids.rds"))
#' labelled_centroids <- get_other_labels(centroids=centroids, tax_level=tax_level, annotation=annotation, MAP=MAP)
#' plot_axes(coordinates=labelled_centroids, tax_level=tax_level, axis1=1, axis2=2, annotation=annotation, MAP=MAP)
plot_axes <- function(coordinates=NULL, centroids, tax_level="ASV", axis1=1, axis2=2, annotation="host",
                      MAP=FALSE, show_plot=FALSE) {
  if(annotation == "host") {
    point_size <- 1
    point_df <- data.frame(ax1=coordinates[coordinates$coord == axis1,]$value,
                          ax2=coordinates[coordinates$coord == axis2,]$value,
                          labels=coordinates[coordinates$coord == axis1,]$labels)
    p <- ggplot() +
      geom_point(data=point_df, aes(x=ax1, y=ax2, color=labels), size=point_size)
    text_df <- data.frame(ax1=centroids[centroids$coord == axis1,]$value,
                          ax2=centroids[centroids$coord == axis2,]$value,
                          labels=centroids[centroids$coord == axis1,]$labels)
    p <- p +
      geom_text(data=text_df, aes(x=ax1, y=ax2, label=labels), color="black", fontface="bold")
    p <- p + theme(legend.position='none')
    img_width <- 4
  } else {
    point_size <- 3
    plot_df <- data.frame(ax1=centroids[centroids$coord == axis1,]$value,
                          ax2=centroids[centroids$coord == axis2,]$value,
                          labels=centroids[centroids$coord == axis1,]$labels)
    p <- ggplot() + geom_point(data=plot_df, aes(x=ax1, y=ax2, color=labels), size=point_size)
    if(annotation == "samplenumber" | annotation == "sampledensity") {
      p <- p + scale_color_gradient(low="blue", high="red")
    }
    img_width <- 4.5
  }
  p <- p + labs(color = annotation)
  
  if(MAP) {
    save_dir <- check_output_dir(c("output","plots",paste0(tax_level,"_MAP")))
  } else {
    save_dir <- check_output_dir(c("output","plots",tax_level))
  }
  plot_save_name <- paste0("Sigma_ordination_",axis1,"x",axis2,"_",annotation,".png")
  
  aspect_ratio.x <- max(centroids[centroids$coord == axis1,]$value) - min(centroids[centroids$coord == axis1,]$value)
  aspect_ratio.y <- max(centroids[centroids$coord == axis2,]$value) - min(centroids[centroids$coord == axis2,]$value)
  img_height <- (aspect_ratio.y/aspect_ratio.x)*img_width
  if(img_height < 2) {
    img_height <- 2
  }
  if(show_plot) {
    show(p)
  }
  ggsave(file.path(save_dir, plot_save_name), plot=p, dpi=150, scale=1.5, width=img_width, height=img_height, units="in")
}

#' Plot first several principle coordinates of ordination
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param which_measure estimated object to visualize, either Lambda or Sigma
#' @param annotation label to assign (e.g. individual)
#' @param MAP use MAP estimate model output instead of full posterior output
#' @param spike_in copy and permute the samples and calculate distances to these too, so as to give upper bound on distances
#' @param show_plot show() plot of first 2 principle coordinates (in addition to rendering first 4 PCoA to files)
#' @return NULL
#' @import ggplot2
#' @export
#' @examples
#' plot_ordination(tax_level="ASV", which_measure="Sigma", annotation="host", MAP=FALSE)
plot_ordination <- function(tax_level="ASV", which_measure="Sigma", annotation="host",
                            MAP=FALSE, spike_in=FALSE, show_plot=FALSE) {
  if(annotation == "host") {
    if(MAP) {
      if(spike_in) {
        coordinates <- read_file(file.path("output","plots",paste0(tax_level,"_MAP"),"Sigma_ordination_spikein.rds"))
        centroids <- read_file(file.path("output","plots",paste0(tax_level,"_MAP"),"Sigma_ordination_centroids_spikein.rds"))
      } else {
        coordinates <- read_file(file.path("output","plots",paste0(tax_level,"_MAP"),"Sigma_ordination.rds"))
        centroids <- read_file(file.path("output","plots",paste0(tax_level,"_MAP"),"Sigma_ordination_centroids.rds"))
      }
    } else {
      coordinates <- read_file(file.path("output","plots",tax_level,"Sigma_ordination.rds"))
      centroids <- read_file(file.path("output","plots",tax_level,"Sigma_ordination_centroids.rds"))
    }
    plot_axes(coordinates=coordinates, centroids=centroids, tax_level=tax_level, axis1=1, axis2=2,
              annotation=annotation, MAP=MAP, show_plot=show_plot)
    if(max(centroids$coord) >= 4) {
      plot_axes(coordinates=coordinates, centroids=centroids, tax_level=tax_level, axis1=3, axis2=4,
                annotation=annotation, MAP=MAP)
    }
  } else {
    if(MAP) {
      # note: ignore "spike-in" directive if we're not labeling by host
      centroids <- read_file(file.path("output","plots",paste0(tax_level,"_MAP"),"Sigma_ordination_centroids.rds"))
    } else {
      centroids <- read_file(file.path("output","plots",tax_level,"Sigma_ordination_centroids.rds"))
    }
    centroids <- get_other_labels(centroids=centroids, tax_level=tax_level, annotation=annotation, MAP=MAP)
    plot_axes(centroids=centroids, tax_level=tax_level, axis1=1, axis2=2, annotation=annotation,
              MAP=MAP, show_plot=show_plot)
    if(max(centroids$coord) >= 4) {
      plot_axes(centroids=centroids, tax_level=tax_level, axis1=3, axis2=4, annotation=annotation,
                MAP=MAP)
    }
    if(max(centroids$coord) >= 6) {
      plot_axes(centroids=centroids, tax_level=tax_level, axis1=5, axis2=6, annotation=annotation,
                MAP=MAP)
    }
    if(max(centroids$coord) >= 8) {
      plot_axes(centroids=centroids, tax_level=tax_level, axis1=7, axis2=8, annotation=annotation,
                MAP=MAP)
    }
  }
}

#' Predict from fitted basset model
#' 
#' @param X observations in days relative to baseline
#' @param fit bassetfit object
#' @param n_samples number of posterior samples to draw
#' @details Predictions will interpolate from beginning to end of observations.
#' @return list of prediction input (days relative to baseline) and output (Eta samples)
#' @import stray
#' @export
#' @examples
#' make_posterior_predictions(X, fit, n_samples=100)
make_posterior_predictions <- function(X, fit) {
  X_predict <- t(1:(max(X)))
  predicted <- predict(fit, X_predict, response="Eta", iter=fit$iter)
  return(list(X_predict=X_predict, Y_predict=predicted))
}

#' Plot MAP covariance matrix for designated host
#' 
#' @param host host short name (e.g. ACA)
#' @param tax_level taxonomic level at which to agglomerate data
#' @param show_plot show() plot of first 2 principle coordinates (in addition to rendering first 4 PCoA to files)
#' @details Output are png files of covariance and correlation matrices.
#' @return NULL
#' @import ggplot2
#' @import driver
#' @export
#' @examples
#' plot_MAP_covariance(host="ZIZ", tax_level="ASV")
plot_MAP_covariance <- function(host, tax_level="ASV", show_plot=FALSE) {
  fit_filename <- file.path("output","model_fits",paste0(tax_level,"_MAP"),paste0(host,"_bassetfit.rds"))
  if(!file.exists(fit_filename)) {
    stop(paste0("MAP model fit for host ",host," at taxonomic level ",tax_level," does not exist!\n"))
  }
  fit_obj <- read_file(fit_filename)
  # plot as covariance
  df <- driver::gather_array(fit_obj$fit$Sigma[,,1], "value", "row", "col")
  p <- ggplot(df, aes(row, col)) +
    geom_tile(aes(fill = value), colour = "white") +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred")
  if(show_plot) {
    show(p)
  }
  save_dir <- check_output_dir(c("output","plots",paste0(tax_level,"_MAP")))
  ggsave(file.path(save_dir,paste0(host,"_MAP_covariance.png")),
         plot=p, scale=1.5, width=5, height=4, units="in", dpi=150)
  # plot as correlation
  df <- driver::gather_array(cov2cor(fit_obj$fit$Sigma[,,1]), "value", "row", "col")
  p <- ggplot(df, aes(row, col)) +
    geom_tile(aes(fill = value), colour = "white") +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred")
  if(show_plot) {
    show(p)
  }
  ggsave(file.path(save_dir,paste0(host,"_MAP_correlation.png")),
         plot=p, scale=1.5, width=5, height=4, units="in", dpi=150)
}

#' Plot posterior predictive intervals of Eta (denoised ALR abundances) over observed time series for designated host
#' 
#' @param host host short name (e.g. ACA)
#' @param tax_level taxonomic level at which to agglomerate data
#' @param predict_coords ALR coordinates to visualize prediction intervals over
#' @param logratio logratio representation to use (e.g. "alr", "ilr", "clr")
#' @param show_plot show() plot of first 2 principle coordinates (in addition to rendering first 4 PCoA to files)
#' @details Output are png files of predicted series intervals.
#' @return NULL
#' @import driver
#' @import dplyr
#' @import ggplot2
#' @import phyloseq
#' @import dplyr
#' @export
#' @examples
#' plot_posterior_predictive(host="ZIB", tax_level="ASV", predict_coords=c(1,2,3))
plot_posterior_predictive <- function(host, tax_level="ASV", predict_coords=NULL, logratio="clr", show_plot=FALSE) {
  fit_filename <- file.path("output","model_fits",tax_level,paste0(host,"_bassetfit.rds"))
  if(!file.exists(fit_filename)) {
    stop(paste0("Model fit for host ",host," at taxonomic level ",tax_level," does not exist!\n"))
  }
  fit_obj <- read_file(fit_filename)
  if(is.null(predict_coords)) {
    # choose random coordinates to predict over
    predict_coords <- sample(1:(fit_obj$fit$D-1))[1:2]
  }

  # a dumb hack for now; these need to be global
  SE_sigma <<- fit_obj$kernelparams$SE_sigma
  SE_rho <<- fit_obj$kernelparams$SE_rho
  PER_sigma <<- fit_obj$kernelparams$PER_sigma
  
  # note: these predictions are in the ALR
  predict_obj <- make_posterior_predictions(fit_obj$X, fit_obj$fit)

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
  
  for(coord in predict_coords) {
    observations <- fit_obj$X
    lr_tidy <- gather_array(lr_ys, "logratio_value", "timepoint", "logratio_coord")
    
    # replace timepoints with observation dates for readability
    map <- data.frame(timepoint=1:length(observations), observation=c(observations))
    lr_tidy <- merge(lr_tidy, map, by="timepoint")
    lr_tidy <- lr_tidy[,!(names(lr_tidy) %in% c("timepoint"))]
    
    no_samples <- dim(Eta)[3]
    posterior_samples <- gather_array(Eta[coord,,], "logratio_value", "observation", "sample_number")
    
    # get quantiles
    post_quantiles <- posterior_samples %>%
      group_by(observation) %>%
      summarise(p2.5 = quantile(logratio_value, prob=0.025),
                p5 = quantile(logratio_value, prob=0.05),
                p10 = quantile(logratio_value, prob=0.1),
                p25 = quantile(logratio_value, prob=0.25),
                p50 = quantile(logratio_value, prob=0.5),
                mean = mean(logratio_value),
                p75 = quantile(logratio_value, prob=0.75),
                p90 = quantile(logratio_value, prob=0.9),
                p95 = quantile(logratio_value, prob=0.95),
                p97.5 = quantile(logratio_value, prob=0.975)) %>%
      ungroup()
    
    p <- ggplot(post_quantiles, aes(x=observation, y=mean)) +
      geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
      geom_ribbon(aes(ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9) +
      geom_line(color="blue") +
      geom_point(data=lr_tidy[lr_tidy$logratio_coord == coord,], aes(x=observation, y=logratio_value), alpha=0.5) +
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle=45)) +
      ylab("LR coord")
    if(show_plot) {
      show(p)
    }
    save_dir <- check_output_dir(c("output","plots",tax_level))
    ggsave(file.path(save_dir,paste0(host,"_posterior_predictive_",coord,".png")),
           plot=p, scale=1.5, width=10, height=2, units="in", dpi=150)
  }
}














