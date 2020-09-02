# Choosing samples for metagenomics
#
# We want to choose samples from among individuals with good "fitness" annotations. We'll brute-force search
# for a combination of individuals that (1) maximizes "dynamics" distance between hosts while (2) making sure
# total number of samples don't exceed some ceiling.

library(ROL)
library(phyloseq)

# --------------------------------------------------------------------------------------------------------------
#   Find hosts with good "fitness" annotations; these are necessarily females!
# --------------------------------------------------------------------------------------------------------------

# raw_data <- readRDS(file.path("input", "ps0.rds"))
data <- load_data(tax_level = "ASV")
metadata <- sample_data(data)
hosts <- unique(sample_data(data)$sname)

outcomes <- read.csv(file.path("input","individual_traits.csv"), header=TRUE)
outcomes <- outcomes[,c("sname", "known_lifespan", "age_first_live_birth", "LRS_livebirths",
                        "LRS_survbirths", "lifetime_rateLiveBirths", "lifetime_rateSurvBirths")]
outcomes <- outcomes[complete.cases(outcomes),]
cat(paste0("There are ",nrow(outcomes)," individuals matching fitness annotation selection criteria!\n"))

hosts <- unique(outcomes$sname)

# get sample count for these individuals
sample_counts <- unlist(as.list(table(metadata$sname)))
sample_counts <- sample_counts[names(sample_counts) %in% hosts]

# eliminate hosts with fewer than 40 samples; we probably won't be able to estimate correlations in these
# with great confidence; only necessary with raw data
# sample_counts <- sample_counts[unname(sample_counts) >= 40]

# --------------------------------------------------------------------------------------------------------------
#   Combinatorially search for the distance maximizing subset of hosts
# --------------------------------------------------------------------------------------------------------------

distances <- readRDS(file.path("output", "Sigma_distance_ASV_MAP.rds"))

include_vec <- distances$host_labels %in% names(sample_counts)
distances$host_labels <- distances$host_labels[include_vec]
d <- distances$distance_mat[include_vec, include_vec]
diag(d) <- NA
# image(d)

# we can define the max number of individuals we'll have to consider as
sample_ceilings <- c(500, 750, 1000)
for(sample_ceiling in sample_ceilings) {
  host_ceiling <- floor(sample_ceiling/min(sample_counts))
  for(m in 2:host_ceiling) {
    start_time <- Sys.time()
    combos <- combn(1:length(sample_counts), m = m)
    combo_stats <- apply(combos, 2, function(x) {
      c(mean(d[x,x], na.rm = TRUE), sum(sample_counts[x]))
    })
    censor_distances <- combo_stats[2,] > sample_ceiling
    combo_stats[1,censor_distances] <- NA
    if(sum(!is.na(combo_stats[1,])) == 0) {
      cat(paste0("No combination of ", m ," hosts that doesn't exceed sample limit!\n"))
      break
    } else {
      max_combo <- which(combo_stats[1,] == max(combo_stats[1,], na.rm = TRUE))
      cat(paste0("Execution: ",round(Sys.time() - start_time, 2)," seconds\n"))
      cat(paste0("Combination of ", m ," hosts that maximizes distances has ",combo_stats[2,max_combo]," / ",sample_ceiling," samples\n"))
      cat("Hosts:", paste(names(sample_counts)[combos[,max_combo]], sep = " "),"\n\n")
    }
  }
}

# --------------------------------------------------------------------------------------------------------------
#   Visualize the selected individuals; obvs. hard-coded
# --------------------------------------------------------------------------------------------------------------

# data <- load_data(tax_level = "ASV")
# metadata <- sample_data(data)
# all_hosts <- unique(metadata$sname)
# for(i in 1:3) {
#   if(i == 1) {
#     selected_hosts <- c("DUN", "LAZ", "OFR", "ONY", "ORI", "VIN") # 500
#   } else if(i == 2) {
#     selected_hosts <- c("DUN", "EAG", "HON", "LAZ", "NOB", "ONY", "ORI", "VIN", "VOT") # 750
#   } else {
#     selected_hosts <- c("DUN", "EAG", "FAX", "HON", "LAZ", "OFR", "ONY", "ORI", "VEL", "VEX", "VIN", "VOT") # 1000
#   }
#   plot_df <- data.frame(sample_date = c(), host = c(), selected = c())
#   group_list <- c()
#   for(host in all_hosts) {
#     # plot_timecourse(data, host = "DUN", gapped = TRUE, legend = FALSE, selected_samples = FALSE, show_plot = FALSE)
#     host_dates <- metadata$collection_date[metadata$sname == host]
#     group_list <- c(group_list, metadata$grp[metadata$sname == host])
#     new_df <- data.frame(sample_date = host_dates, host = host, selected = (host %in% selected_hosts))
#     plot_df <- rbind(plot_df, new_df)
#   }
#   plot_df$sample_date <- as.Date(plot_df$sample_date)
#   p <- ggplot(plot_df, aes(x = sample_date, y = host, color = selected)) +
#   geom_point()
#   ggsave(paste0("set_",i,".png"), p, units = "in", dpi = 150, height = 6, width = 6)
#   # get unique groups represented by these hosts
#   print(unique(group_list))
# }

# --------------------------------------------------------------------------------------------------------------
#   Greedily find distance maximizers
# --------------------------------------------------------------------------------------------------------------

# distances <- readRDS(file.path("output", "Sigma_distance_ASV_MAP.rds"))
 
# # we'll use a greedy solution; the combinatorics here are awful
# # alternatively, we could think about doing a constrained optimization (and I did) to optimize a vector
# # of indicators that maximizes distances between hosts subject to a constraint on total sample number
# # but it's... non trivial :)
# 
# select_hosts <- function(distances, sample_counts, sample_ceiling) {
#   # the jist: find the most distant individuals, then find the individual the farthest from those, then
#   # find the individual farthest from those... etc.
#   # exclude distances for hosts not selected on the basis of fitness annotations
#   exclude_vec <- !(distances$host_labels %in% hosts)
#   d <- distances$distance_mat
#   # d[upper.tri(d, diag = T)] <- NA
#   diag(d) <- NA
#   d[exclude_vec,] <- NA
#   d[,exclude_vec] <- NA
#   # image(d)
#   
#   total_samples <- 0
#   hosts_selected <- numeric(length(sample_counts))
#   names(hosts_selected) <- names(sample_counts)
#   max_pair <- which(d == max(d, na.rm = TRUE), arr.ind = TRUE)[1,]
#   h1 <- distances$host_labels[max_pair[["row"]]]
#   h2 <- distances$host_labels[max_pair[["col"]]]
#   d[max_pair[["row"]], max_pair[["col"]]] <- NA
#   d[max_pair[["col"]], max_pair[["row"]]] <- NA
#   hosts_selected[[h1]] <- 1
#   hosts_selected[[h2]] <- 1
#   total_samples <- sum(sample_counts[which(hosts_selected == 1)])
#   repeat {
#     # find the next, jointly farthest individual from the ones we've already selected
#     search_idx <- which(distances$host_labels %in% names(hosts_selected[hosts_selected == 1]))
#     
#     # if(sum(!is.na(d[search_idx,])) == 0) {
#     #   # we couldn't find anything that maximizes the distance to the individuals already selected
#     #   # that *doesn't* exceed the sample number limit
#     #   # expand the search to all
#     #   search_idx <- rep(TRUE, nrow(d))
#     # }
#     
#     joint_rows <- apply(d[,search_idx], 1, mean)
#     max_host <- which(joint_rows == max(joint_rows, na.rm = TRUE))
#     d[max_host,search_idx] <- NA
#     d[search_idx,max_host] <- NA
#     new_h <- distances$host_labels[max_host]
#     new_hosts_selected <- hosts_selected
#     new_hosts_selected[[new_h]] <- 1
#     new_total_samples <- sum(sample_counts[which(new_hosts_selected == 1)])
#     if(new_total_samples > sample_ceiling) {
#       margin <- sample_ceiling - total_samples
#       min_remaining <- min(sample_counts[names(sample_counts) %in% names(hosts_selected[hosts_selected == 0])])
#       if(margin < min_remaining) {
#         break
#       } # else continue
#     } else {
#       hosts_selected <- new_hosts_selected
#       total_samples <- new_total_samples
#     }
#   }
#   
#   return(list(hosts = names(hosts_selected[hosts_selected == 1]), total_samples = total_samples))
# }
# 
# res <- select_hosts(distances, sample_counts, sample_ceiling = 500)
# res <- select_hosts(distances, sample_counts, sample_ceiling = 750)
# res <- select_hosts(distances, sample_counts, sample_ceiling = 1000)
# 
# # we could think of improving this by reaching for the next most distant


