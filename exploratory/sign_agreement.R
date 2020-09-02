library(ROL)

# get percent of shared sign in correlation across individuals
Sigmas <- load_full_posteriors(tax_level = "ASV", logratio = "clr")
k <- dim(Sigmas[[1]])[3]
m <- length(Sigmas)
agree_percent <- numeric(k)
for(i in 1:k) {
  cat("Posterior sample",i,"\n")
  correlations.1 <- lapply(Sigmas, function(x) {
    temp <- cov2cor(x[,,i])
    temp[upper.tri(temp, diag = FALSE)]
  })
  correlations.2 <- matrix(unlist(correlations.1), ncol = length(correlations.1[[1]]), byrow = TRUE)
  signs <- sign(correlations.2)
  temp <- apply(signs, 2, function(x) {
    max(table(x))/m
  })
  agree_percent[i] <- mean(temp)
}

round(agree_percent, 2)
