#' Calculate Riemannian distances between all posterior samples
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param which_measure estimated object to embed, either Lambda or Sigma
#' @param MAP use MAP estimate model output instead of full posterior output
#' @details Distance matrix between posterior samples saved to designated output directory. Saves
#' a list containing distance matrix and labels of each row or column (by host).
#' @return NULL
#' @import driver
#' @export
#' @examples
#' calc_posterior_distances(tax_level="genus", which_measure="Sigma", MAP=FALSE)
calc_posterior_distances <- function(tax_level="genus", which_measure="Sigma", MAP=FALSE) {
  # grab all fitted models
  model_list <- get_fitted_model_list(tax_level=tax_level, MAP=MAP)
  P <- model_list$D - 1 # ALR
  n_samples <- model_list$n_samples
  n_hosts <- length(model_list$hosts)
  # initialize samples matrix
  if(which_measure == "Sigma") {
    all_samples <- matrix(NA, P, P*n_samples*n_hosts)
  } else {
    # we'll use per-sample average to manage individuals having different N so each individual's
    # posterior will be summarized as one mean vector
    all_samples <- matrix(NA, n_hosts*n_samples, P)
  }
  # insert samples (column-wise) into samples matrix
  host_labels <- c()
  for(i in 1:n_hosts) {
    fit <- fix_MAP_dims(readRDS(model_list$model_list[i])$fit)
    # convert to ILR; this can be removed
    V <- driver::create_default_ilr_base(ncategories(fit))
    fit.ilr <- to_ilr(fit, V)
    Lambda <- fit.ilr$Lambda
    Sigma <- fit.ilr$Sigma
    if(which_measure == "Sigma") {
      Sigma <- Sigma[,,1:n_samples]
      all_samples[,((i-1)*P*n_samples+1):(i*P*n_samples)] <- Sigma
      host_labels <- c(host_labels, rep(model_list$hosts[i], n_samples))
    } else {
      collLambda <- t(apply(Lambda, 3, function(X) { apply(X, 1, mean) })) # n_samples x P
      all_samples[((i-1)*n_samples+1):(i*n_samples),] <- collLambda
      host_labels <- c(host_labels, rep(model_list$hosts[i], n_samples))
    }
  }
  
  save_dir <- check_output_dir(c("output"))
  if(MAP) {
    dist_filename <- file.path(save_dir,paste0(which_measure,"_distance_",tax_level,"_MAP.rds"))
  } else {
    dist_filename <- file.path(save_dir,paste0(which_measure,"_distance_",tax_level,".rds"))
  }
  if(which_measure == "Sigma") {
    distance_mat <- Riemann_dist_samples(all_samples, n_hosts, n_samples)
  } else {
    # use Euclidean distance
    distance_mat <- as.matrix(dist(all_samples))
  }
  saveRDS(list(host_labels=host_labels, distance_mat=distance_mat), file=dist_filename)
}

#' Embeds posterior samples using MDS and a pre-calculated distance matrix
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param which_measure estimated object to embed, either Lambda or Sigma
#' @param MAP use MAP estimate model output instead of full posterior output
#' @details Distance matrix between posterior samples must be present in designated output directory
#' @return NULL
#' @import ggplot2
#' @import tidyr
#' @export
#' @examples
#' embed_posteriors(tax_level="genus", which_measure="Sigma", MAP=FALSE)
embed_posteriors <- function(tax_level="genus", which_measure="Sigma", MAP=FALSE) {
  if(MAP) {
    dist_filename <- file.path("output",paste0(which_measure,"_distance_",tax_level,"_MAP.rds"))
  } else {
    dist_filename <- file.path("output",paste0(which_measure,"_distance_",tax_level,".rds"))
  }
  if(!file.exists(dist_filename)) {
    stop(paste0("Error: unable to locate distance matrix over posterior samples at level ",tax_level,"!\n"))
  }
  dist_obj <- readRDS(dist_filename)
  host_labels <- dist_obj$host_labels
  distance_mat <- dist_obj$distance_mat
  
  fit <- cmdscale(distance_mat, eig=TRUE, k=100)
  # I believe in this case the magnitude of the eigenvalues is proportional to the variance
  # explained (technically it differs by a factor of n-1 [the DOF] I think)
  eig_tot <- sum(abs(fit$eig))
  cat(paste0("Eigenvalue #1: ",round(fit$eig[1],2)," (% variance: ",round(abs(fit$eig[1])/eig_tot,2),")\n"))
  cat(paste0("Eigenvalue #2: ",round(fit$eig[2],2)," (% variance: ",round(abs(fit$eig[2])/eig_tot,2),")\n"))
  cat(paste0("Eigenvalue #3: ",round(fit$eig[3],2)," (% variance: ",round(abs(fit$eig[3])/eig_tot,2),")\n"))
  cat(paste0("Eigenvalue #4: ",round(fit$eig[4],2)," (% variance: ",round(abs(fit$eig[4])/eig_tot,2),")\n"))
  cat(paste0("Eigenvalue #5: ",round(fit$eig[5],2)," (% variance: ",round(abs(fit$eig[5])/eig_tot,2),")\n"))
  cat(paste0("Eigenvalue #6: ",round(fit$eig[6],2)," (% variance: ",round(abs(fit$eig[6])/eig_tot,2),")\n"))
  
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
    saveRDS(df, file.path("output","plots",paste0(tax_level,"_MAP"),paste0(which_measure,"_ordination.rds")))
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
    saveRDS(df_centroids, file.path("output","plots",paste0(tax_level,"_MAP"),paste0(which_measure,"_ordination_centroids.rds")))
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
#' labelled_centoids <- get_other_labels(centroids, tax_level="genus", annotation="group", MAP=FALSE)
get_other_labels <- function(centroids, tax_level="genus", annotation="group", MAP=FALSE) {
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
    for(host in host) {
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
#' @details Coordinates parameter is only necessary where we're plotting full host posteriors (i.e. where
#' annotation="individual").
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
plot_axes <- function(coordinates=NULL, centroids, tax_level="genus", axis1=1, axis2=2, annotation="individual", MAP=FALSE) {
  if(annotation == "individual") {
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
  ggsave(file.path(save_dir, plot_save_name), plot=p, dpi=150, scale=1.5, width=img_width, height=img_height, units="in")
}

#' Plot a pair of principal coordinates
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param which_measure estimated object to visualize, either Lambda or Sigma
#' @param axis2 PCoA coordinate to display on y-axis
#' @param annotation label to assign (e.g. individual)
#' @param MAP use MAP estimate model output instead of full posterior output
#' @return NULL
#' @import ggplot2
#' @export
#' @examples
#' plot_ordination(tax_level="genus", axis1=1, axis2=2, annotation="individual", MAP=FALSE)
plot_ordination <- function(tax_level="genus", which_measure="Sigma", annotation="individual", MAP=FALSE) {
  if(annotation == "individual") {
    if(MAP) {
      coordinates <- readRDS(file.path("output","plots",paste0(tax_level,"_MAP"),"Sigma_ordination.rds"))
      centroids <- readRDS(file.path("output","plots",paste0(tax_level,"_MAP"),"Sigma_ordination_centroids.rds"))
    } else {
      coordinates <- readRDS(file.path("output","plots",tax_level,"Sigma_ordination.rds"))
      centroids <- readRDS(file.path("output","plots",tax_level,"Sigma_ordination_centroids.rds"))
    }
    plot_axes(coordinates=coordinates, centroids=centroids, tax_level=tax_level, axis1=1, axis2=2, annotation=annotation, MAP=FALSE)
    plot_axes(coordinates=coordinates, centroids=centroids, tax_level=tax_level, axis1=3, axis2=4, annotation=annotation, MAP=FALSE)
  } else {
    if(MAP) {
      centroids <- readRDS(file.path("output","plots",paste0(tax_level,"_MAP"),"Sigma_ordination_centroids.rds"))
    } else {
      centroids <- readRDS(file.path("output","plots",tax_level,"Sigma_ordination_centroids.rds"))
    }
    centroids <- get_other_labels(centroids=centroids, tax_level=tax_level, annotation=annotation, MAP=MAP)
    plot_axes(centroids=centroids, tax_level=tax_level, axis1=1, axis2=2, annotation=annotation, MAP=FALSE)
    plot_axes(centroids=centroids, tax_level=tax_level, axis1=3, axis2=4, annotation=annotation, MAP=FALSE)
  }
}




