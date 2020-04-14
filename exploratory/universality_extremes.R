library(ROL)
library(driver)
library(ggplot2)

heatmap <- plot_interaction_heatmap(tax_level="ASV", logratio="clr", Sigmas=NULL,
                                         taxon_idx=NULL, show_plot=FALSE, return_matrix=TRUE)

p <- nrow(heatmap)
q <- ncol(heatmap)
sign_matrix <- matrix(sapply(c(heatmap), sign), p, q)

visualize <- FALSE

if(visualize) {
  df <- gather_array(sign_matrix, "sign", "host", "interaction")
  p <- ggplot(df, aes(interaction, host)) +
        geom_tile(aes(fill = sign)) +
        scale_fill_gradient2(low = "darkblue", high = "darkred")
  ggsave("test.png", p, units="in", dpi=100, height=5, width=15)

  # permute and visualize
  permuted_matrix <- t(apply(sign_matrix, 1, function(x) x[sample(1:q)]))

  df <- gather_array(permuted_matrix, "sign", "host", "interaction")
  p <- ggplot(df, aes(interaction, host)) +
        geom_tile(aes(fill = sign)) +
        scale_fill_gradient2(low = "darkblue", high = "darkred")
  ggsave("test2.png", p, units="in", dpi=100, height=5, width=15)
}

# score

# mat is a sign matrix
# let the score be the number of same-sign interactions within each host, summed
score <- function(mat) {
  score <- 0
  for(i in 1:q) {
    if(sum(mat[,i]) < 0) {
      # negative correlations dominate
      score <- score + sum((mat[,i]*(-1) + 1) / 2)
    } else {
      # positive correlations dominate
      score <- score + sum((mat[,i] + 1) / 2)
    }
  }
  return(score)
}

iterations <- 100000

p_scores <- c()
for(j in 1:iterations) {
  if(iterations %% 1000 == 0) {
    cat("Iteration",j,"\n")
  }
  permuted_matrix <- t(apply(sign_matrix, 1, function(x) x[sample(1:q)]))
  p_scores <- c(p_scores, score(permuted_matrix))
}

true_score <- score(sign_matrix)

saveRDS(list(p_scores=p_scores, true_score=true_score), "p_scores.rds")

