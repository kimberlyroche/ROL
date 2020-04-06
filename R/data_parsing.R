#' Load file with existence check
#' 
#' @param filename filename (and path)
#' @return RDS file
#' @export
#' @examples
#' input_data <- read_file("input/filtered_genus_5_20.rds")
read_file <- function(filename) {
  if(!file.exists(filename)) {
    stop(paste0("No such file ",filename,"!\n"))
  }
  return(readRDS(filename))
}

#' Load ABRP 16S data and filter out low abundance taxa
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param host_sample_min minimum sample number for host inclusion in the filtered data set
#' @param count_threshold minimum count at which a taxon must be observed
#' @param sample_threshold minimum proportion of samples within each host at which a taxon must be observed at the level specified by count_threshold
#' @details Together count_threshold and sample_threshold specify a minimum representation for a taxon.
#' Taxa below these observed thresholds will be grouped together into a category "Other."
#' @return phyloseq object
#' @import phyloseq
#' @export
#' @examples
#' data <- load_data(tax_level="ASV", host_sample_min=75, count_threshold=5, sample_threshold=0.2)
load_data <- function(tax_level="ASV", host_sample_min=75, count_threshold=5, sample_threshold=0.2) {
  if(is.null(tax_level)) {
    tax_level <- "ASV"
  }
  if(tax_level == "ASV") {
    # check to see if a pre-filtered version of this is available
    filename <- file.path("input", paste0("filtered_",tax_level,"_",count_threshold,"_",round(sample_threshold*100),".rds"))
    if(file.exists(filename)) {
      return(readRDS(filename))
    }
    # continue with no agglomeration
    filename <- file.path("input","ps0.rds")
    agglomerated_data <- readRDS(filename)
  } else {
    # agglomerate
    agglomerated_data <- agglomerate_data(tax_level=tax_level)
  }
  # subset to hosts with minimum sample number
  sname_occurrences <- sample_data(agglomerated_data)$sname
  # global assign is a hack seemingly necessary for this phyloseq::subset_samples function call
  hosts_over_threshold <<- names(which(table(sname_occurrences) >= host_sample_min))
  subsetted_data <- subset_samples(agglomerated_data, sname %in% hosts_over_threshold)
  # filter
  filtered_data <- filter_data(subsetted_data, tax_level=tax_level, count_threshold=count_threshold, sample_threshold=sample_threshold)
  return(filtered_data)
}

#' Agglomerate ABRP 16S data to desired taxonomic level
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @return phyloseq object
#' @import phyloseq
#' @importFrom phyloseq tax_glom
#' @export
#' @examples
#' data <- agglomerate_data(tax_level="ASV")
agglomerate_data <- function(tax_level="ASV") {
  if(tax_level == "ASV") {
    return(data)
  }
  out_filename <- file.path("input",paste0("glom_data_",tax_level,".rds"))
  if(file.exists(out_filename)) {
    data <- readRDS(out_filename)
    return(data)
  }
  in_filename <- file.path("input","ps0.rds")
  if(file.exists(in_filename)) {
    data <- readRDS(in_filename)
  } else {
    stop(paste0("Input data file ",in_filename," does not exist!\n"))
  }
  agglomerated_data <- tax_glom(data, taxrank=tax_level, NArm=FALSE)
  saveRDS(agglomerated_data, file=out_filename)
  return(agglomerated_data)
}

#' Filter out low abundance taxa
#' 
#' @param data a phyloseq object
#' @param tax_level taxonomic level at which data have been agglomerated
#' @param count_threshold minimum count at which a taxon must be observed
#' @param sample_threshold minimum proportion of samples within each host at which a taxon must be observed at the level specified by count_threshold
#' @details Together count_threshold and sample_threshold specify a minimum representation for a taxon.
#' Taxa below these observed thresholds will be grouped together into a category "Collapsed".
#' @return phyloseq object
#' @import phyloseq
#' @export
#' @examples
#' filtered_data <- filter_data(data, count_threshold=5, sample_threshold=0.2)
filter_data <- function(data, tax_level=NULL, count_threshold=5, sample_threshold=0.2) {
  if(is.null(tax_level)) {
    stop("Error: omitted taxonomic level label in filter_data()!\n")
  }
  filename <- file.path("input", paste0("filtered_",tax_level,"_",count_threshold,"_",round(sample_threshold*100),".rds"))
  if(!file.exists(filename)) {
    snames <- unique(sample_data(data)$sname)
    count_table <- otu_table(data)@.Data # samples x taxa
    total_counts <- sum(count_table)
    # iterate each individual, failing taxa that fall below the specified thresholds in any individuals
    keep_indices <- rep(TRUE, ntaxa(data))
    for(host in snames) {
      host <<- host
      subset_data <- subset_samples(data, sname == host)
      subset_count_table <- otu_table(subset_data)@.Data
      # these are the indices to remove!
      host_keep_indices <- as.vector(apply(subset_count_table, 2, function(x) sum(x >= count_threshold)/phyloseq::nsamples(subset_data) >= sample_threshold))
      keep_indices <- keep_indices & host_keep_indices
    }
    collapse_indices <- !keep_indices
    # collapse Mitochondria, Chloroplasts, and anything outside the bacterial domain
    # we'll retain these counts, as they inform sampling depth and so uncertainty, but we're not interested in their interactions
    tt <- tax_table(data)@.Data
    collapse_indices[which(tt[,colnames(tt) == "family"] == "Mitochondria")] <- TRUE
    collapse_indices[which(tt[,colnames(tt) == "order"] == "Chloroplast")] <- TRUE
    # exclude Archaea or not?
    # non_bacterial_taxa <- collapse_indices[which(tt[,colnames(tt) == "domain"] != "Bacteria")]
    # if(length(non_bacterial_taxa) > 0) {
    #   collapse_indices[non_bacterial_taxa] <- TRUE
    # }
    # perform the collapsing via phyloseq
    merged_data <- merge_taxa(data, which(collapse_indices == TRUE), 1)
    retained_counts <- count_table[,!collapse_indices]
    cat(paste0("Retaining ",sum(collapse_indices == FALSE)," / ",ntaxa(data)," taxa\n"))
    cat("Collapsed ",round((total_counts-sum(retained_counts))/total_counts, 3)*100,"% of total counts\n")
    # calculate proportion zeros
    cat("Percent zeros is ",round(sum(retained_counts == 0)/(nrow(retained_counts)*ncol(retained_counts)), 3)*100,"%\n")
    # collapse all into first index; we'll label this domain : species "Collapsed"
    collapsed_into_idx <- which(collapse_indices == TRUE)[1]
    cat(paste0("Other category is index ",collapsed_into_idx,"\n"))
    tax_table(merged_data)@.Data[collapsed_into_idx,] <- rep("Collapsed", 7)
    # save results
    saveRDS(merged_data, file=filename)
    return(merged_data)
  } else {
    return(read_file(filename))
  }
}

#' Create output directory/ies that may not already exist
#' 
#' @param path_to_dir a character vector of subdirectories to create (if not already in existence)
#' @return the path as a joined string
#' @export
#' @examples
#' check_output_dir(list("output", "plots"))
check_output_dir <- function(path_to_dir) {
  for(i in 1:length(path_to_dir)) {
    dir.create(do.call(file.path, as.list(path_to_dir[1:i])), showWarnings=FALSE)
  }
  return(do.call(file.path, as.list(path_to_dir)))
}

#' Get taxonomy, reordered to reflect altered ALR reference
#' 
#' @param data a phyloseq object
#' @param alr_ref index of ALR reference taxon
#' @return data.frame of re-ordered taxonomy
#' @import phyloseq
#' @export
#' @examples
#' tax <- get_taxonomy(data, alr_ref=67)
get_taxonomy <- function(data, alr_ref) {
  if(is.null(alr_ref)) {
    stop(paste0("Missing ALR reference index!\n"))
  }
  if(alr_ref < 1 | alr_ref > ntaxa(data)) {
    stop("Invalid ALR reference index!\n")
  }
  tax <- tax_table(data)@.Data
  tax <- tax[c(setdiff(1:nrow(tax),alr_ref),alr_ref),]
  rownames(tax) <- NULL
  return(tax)
}

#' Get highest taxonomic level to which a microbe is resolved (i.e. not NA)
#' 
#' @param taxonomy a named list of taxonomic assignment for a designated sequence variant
#' @param deepest_tax_level the finest resolution taxonomy level to return (if possible)
#' @return readable label of the form "phylum Tenericutes"
#' @export
#' @examples
#' data <- load_data(tax_level="ASV")
#' alr_ref <- formalize_parameters(data)$alr_ref
#' taxonomy <- get_taxonomy(data, alr_ref)
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

#' Get concise taxonomic labels for microbes
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param logratio logratio representation to use (e.g. "alr", "clr", "none")
#' @param other_category "Other" taxonomic category of collapsed low abundance taxa (if available)
#' @details If other_category is provided, these taxa will be rendered as *LR(Other).
#' @return readable list of taxonomic labels (e.g. of the form "CLR(phylum Tenericutes)")
#' @export
#' @examples
#' tax_labels <- assign_concise_taxonomy(tax_level="ASV", logratio="alr", other_category=0)
assign_concise_taxonomy <- function(tax_level="ASV", logratio="alr", other_category=0) {
  if(logratio != "alr" & logratio != "clr" & logratio != "none") {
    stop(paste0("Only logratio representations ALR and CLR allowed!\n"))
  }
  data <- read_file(file.path("input",paste0("filtered_",tax_level,"_5_20.rds")))
  alr_ref <- formalize_parameters(data)$alr_ref
  taxonomy <- get_taxonomy(data, alr_ref)
  level_number <- which(colnames(taxonomy) == tax_level)
  if(logratio == "alr") {
    labels <- character(nrow(taxonomy)-1)
    alr.ref <- get_deepest_assignment(taxonomy[nrow(taxonomy),])
  } else {
    labels <- character(nrow(taxonomy))
  }
  for(i in 1:length(labels)) {
    if(i == other_category) {
      if(logratio == "alr") {
        labels[i] <- paste0("ALR(Other / ",alr.ref,")")
      } else if(logratio == "clr") {
        labels[i] <- paste0("CLR(Other)")
      } else {
        labels[i] <- paste0("Other")
      }
    } else {
      deep_label <- get_deepest_assignment(taxonomy[i,])
      if(logratio == "alr") {
        if(is.null(deep_label)) {
          labels[i] <- paste0("ALR(Unresolved/",alr.ref,")")
        } else {
          labels[i] <- paste0("ALR(",get_deepest_assignment(taxonomy[i,]),"/",alr.ref,")")
        }
      } else if(logratio == "clr") {
        if(is.null(deep_label)) {
          labels[i] <- paste0("CLR(Unresolved)")
        } else {
          labels[i] <- paste0("CLR(",get_deepest_assignment(taxonomy[i,]),")")
        }
      } else {
        if(is.null(deep_label)) {
          labels[i] <- "Unresolved"
        } else {
          labels[i] <- get_deepest_assignment(taxonomy[i,])
        }
      }
    }
  }
  return(labels)
}
