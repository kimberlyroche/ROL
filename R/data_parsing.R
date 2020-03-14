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
#' @param replicates include replicates
#' @param count_threshold minimum count at which a taxon must be observed
#' @param sample_threshold minimum proportion of samples within each host at which a taxon must be observed at the level specified by count_threshold
#' @details Together count_threshold and sample_threshold specify a minimum representation for a taxon.
#' Taxa below these observed thresholds will be grouped together into a category "Other."
#' @return phyloseq object
#' @import phyloseq
#' @export
#' @examples
#' data <- load_data(tax_level="genus", replicates=TRUE, count_threshold=5, sample_threshold=0.2)
load_data <- function(tax_level="genus", replicates=TRUE, count_threshold=5, sample_threshold=0.2) {
  if(is.null(tax_level)) {
    tax_level <- "ASV"
  }
  if(tax_level == "ASV") {
    # no agglomeration
    if(replicates) {
      filename <- file.path("input","data_w_rep.rds")
    } else {
      filename <- file.path("input","data_wo_rep.rds")
    }
  } else {
    # agglomerate
    agglomerated_data <- agglomerate_data(tax_level=tax_level, replicates=replicates)
  }
  # subset to hosts with minimum sample number
  # global assign is a hack seemingly necessary for this phyloseq::subset_samples function call
  hosts_over_40 <<- c("DUI", "ECH", "LOG", "VET", "DUX", "LEB", "ACA", "OPH", "THR", "VAI", "VIG", "VOG", "DAS",
                     "CAI", "COB", "PEB", "OXY", "WRI", "NAP", "SEB", "COO", "LAD", "LOB", "WAD", "GAB", "LIW",
                     "VIN", "TAL", "VEX", "VEI", "ALE", "MBE", "WHE", "WYN", "LOL", "HOL", "NOB", "VOT", "LYE",
                     "HON", "DAG", "DUN", "OTI", "LUI", "OFR", "LAZ", "ONY", "VEL", "ELV", "FAX", "ORI", "EAG",
                     "ODE", "NIK", "VAP", "WIP", "LOU", "NOO", "EVA", "EXO", "KOR", "NAR", "VOW", "HYM", "PAI",
                     "LAS", "VIO", "WEA", "DOU", "LIZ", "WAS", "ZIB", "QUA", "WEN", "WOB", "WOL", "HOK", "LAV",
                     "OBI", "POK", "SOR", "KOL", "ISR", "OMO", "SCE", "AFR", "MON", "NIN", "VEB", "ADD", "VOY",
                     "DRO", "LOC", "OJU", "OST", "DUB", "LEI", "VAA", "GAN", "HUM", "LUN", "VIV", "BUC", "LAN",
                     "LOX", "HAS", "SNA", "WUA", "YAI", "EGO", "ABB", "CRU", "LOF", "WAB", "ZIZ", "COD", "LEX",
                     "RAJ", "KIW", "LAO", "LIB", "NJU", "OBR", "OCE", "POW", "IAG", "MLO", "GYP", "LIT", "OPA",
                     "COT", "DIP", "LAW", "RHO", "VOR", "AMA", "AYU", "DUR", "FLA", "OAS", "VIB", "CAB", "CHE",
                     "HAV", "LUP", "MIC", "YOB", "PIT", "YAN", "LOZ", "TOG", "BEA", "DUD", "GOM", "HIB", "WAG",
                     "ETO", "KEL", "NUT", "WES", "IDI", "ISO", "PRU", "YOG", "ZAI")
  subsetted_data <- subset_samples(agglomerated_data, sname %in% hosts_over_40)
  # filter
  filtered_data <- filter_data(subsetted_data, tax_level=tax_level, count_threshold=count_threshold, sample_threshold=sample_threshold)
  return(filtered_data)
}

#' Agglomerate ABRP 16S data to desired taxonomic level
#' 
#' @param tax_level taxonomic level at which to agglomerate data
#' @param replicates include replicates
#' @return phyloseq object
#' @import phyloseq
#' @importFrom phyloseq tax_glom
#' @export
#' @examples
#' data <- agglomerate_data(tax_level="genus", replicates=TRUE)
agglomerate_data <- function(tax_level="genus", replicates=TRUE) {
  if(tax_level == "ASV") {
    return(data)
  }
  if(replicates) {
    out_filename <- file.path("input",paste0("glom_data_",tax_level,"_reps.rds"))
  } else {
    out_filename <- file.path("input",paste0("glom_data_",tax_level,".rds"))
  }
  if(file.exists(out_filename)) {
    data <- readRDS(out_filename)
    return(data)
  }
  if(replicates) {
    in_filename <- file.path("input","data_w_rep.rds")
  } else {
    in_filename <- file.path("input","data_wo_rep.rds")
  }
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
    retained_counts <- sum(count_table[,!collapse_indices])
    cat(paste0("Retaining ",sum(collapse_indices == FALSE)," / ",ntaxa(data)," taxa\n"))
    cat("Collapsed ",round((total_counts-retained_counts)/total_counts, 3)*100,"% of total counts\n")
    # calculate proportion zeros
    cat("Percent zeros is ",round(sum(count_table == 0)/(nrow(count_table)*ncol(count_table)), 3)*100,"%\n")
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
