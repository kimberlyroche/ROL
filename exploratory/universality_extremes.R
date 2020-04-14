library(ROL)
library(driver)
library(ggplot2)

heatmap <- plot_interaction_heatmap(tax_level="ASV", logratio="clr", Sigmas=NULL,
                                         taxon_idx=NULL, show_plot=FALSE, return_matrix=TRUE)

discretize <- FALSE

plot_heatmap <- function(mat, filename) {
  df <- gather_array(mat, "sign", "host", "interaction")
  p <- ggplot(df, aes(interaction, host)) +
        geom_tile(aes(fill = sign)) +
        scale_fill_gradient2(low = "darkblue", high = "darkred")
  ggsave(filename, p, units="in", dpi=100, height=5, width=15)  
}

score <- function(mat, discrete=FALSE) {
  # let the score be the number of same-sign interactions within each host, summed
  score_val <- 0
  if(discrete) {
    for(i in 1:ncol(mat)) {
      if(sum(mat[,i]) < 0) {
        # negative correlations dominate
        score_val <- score_val + sum((mat[,i]*(-1) + 1) / 2)
      } else {
        # positive correlations dominate
        score_val <- score_val + sum((mat[,i] + 1) / 2)
      }
    }
  } else {
    score_val <- sum(apply(mat, 2, sd))
  }
  return(score_val)
}

permute_and_score <- function(mat, iterations=1000, discrete=FALSE) {
  # permute and test "significance"
  print_it <- round(iterations/10)
  p_scores <- c()
  for(j in 1:iterations) {
    if(j %% print_it == 0) {
      cat("Iteration",j,"/",iterations,"\n")
    }
    permuted_matrix <- t(apply(mat, 1, function(x) x[sample(1:length(x))]))
    p_scores <- c(p_scores, score(permuted_matrix, discrete=discrete))
  }
  result <- list(p_scores=p_scores, true_score=score(mat, discrete=discrete))
  saveRDS(result, "p_scores.rds")
  return(result)
}

if(discretize) {
  sign_matrix <- matrix(sapply(c(heatmap), sign), p, q)
  score_val <- score(sign_matrix, discrete=TRUE)

  # permute and visualize
  permuted_matrix <- t(apply(sign_matrix, 1, function(x) x[sample(1:length(x))]))
  plot_heatmap(permuted_matrix, "permute_discrete.png")

  # permute and test "significance"
  #p_scores <- permute_and_score(sign_matrix, iterations=100, discrete=TRUE)
} else {
  score_val <- score(heatmap, discrete=FALSE)

  # permute and visualize
  permuted_matrix <- t(apply(heatmap, 1, function(x) x[sample(1:length(x))]))

  plot_heatmap(permuted_matrix, "permute_continuous.png")

  # permute and test "significance"
  #p_scores <- permute_and_score(heatmap, iterations=100, discrete=TRUE)
}



