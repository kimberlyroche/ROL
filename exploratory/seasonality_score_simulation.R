library(ggplot2)
library(lomb)

downsample_series <- function(y1, y2, downsample_to_percent = 0.1) {
  n <- length(y1)
  sample_idx <- sort(sample(1:n)[1:round(n*downsample_to_percent)])
  y1.down <- rep(NA, n)
  y1.down[sample_idx] <- y1[sample_idx]
  y2.down <- rep(NA, n)
  y2.down[sample_idx] <- y2[sample_idx]
  list(y1 = y1.down, y2 = y2.down)
}

plot_series <- function(x, y1, y2, labels = NULL) {
  n <- length(x)
  if(is.null(labels)) {
    labels <- as.factor(c(rep(1, n), rep(2, n)))
  }
  df <- data.frame(day = x, abundance = c(y1, y2), label = labels)
  p <- ggplot(df, aes(x = day, y = abundance, color = label)) +
    geom_point()
  p
}

# labels are seasonal "wet"/"dry" labels (etc.)
score <- function(y1, y2, labels) {
  mean.y1 <- mean(y1)
  mean.y2 <- mean(y2)
  binary.y1 <- sapply(y1, function(yy) {
    if(yy >= mean.y1) {
      1
    } else {
      -1
    }
  })
  binary.y2 <- sapply(y2, function(yy) {
    if(yy >= mean.y2) {
      1
    } else {
      -1
    }
  })
  n <- length(labels)
  series_consensus <- sapply(1:n, function(yy) {
    if(binary.y1[yy] == binary.y2[yy]) {
      if(binary.y1[yy] == 1) {
        1
      } else {
        -1
      }
    } else {
      0
    }
  })
  season_disagreement <- 0
  wet_season <- which(labels == "wet")
  dry_season <- which(labels == "dry")
  wet_total <- sum(series_consensus[wet_season])
  dry_total <- sum(series_consensus[dry_season])
  max_difference <- length(wet_season) + length(dry_season)
  if(sign(wet_total) == sign(dry_total)) {
    observed_difference <- abs(abs(wet_total) - abs(dry_total))
  } else {
    observed_difference <- abs(abs(wet_total) + abs(dry_total))
  }
  observed_difference / max_difference
}

## -----------------------------------------------------------------------------------------------------------------------------
##   Correlated series y, y.corr, anti-correlated series y.acorr, and uncorrelated series y.uncorr
## -----------------------------------------------------------------------------------------------------------------------------

x <- seq(from = 0, to = 10*pi, length.out = 1000)
y <- sin(x) + rnorm(1000)*0.5
y.corr <- sin(x) + rnorm(1000)*0.5
y.acorr <- sin(-x) + rnorm(1000)*0.5
y.uncorr <- rnorm(1000)*0.5
labels <- as.factor(sin(x) > 0)
levels(labels) <- c("wet", "dry")
plot_series(x, y, y.corr)
plot_series(x, y, y.acorr)
plot_series(x, y, y.uncorr)
plot_series(x, y, y.corr, labels = labels)

## -----------------------------------------------------------------------------------------------------------------------------
##   "Score" seasonal co-occurrence
## -----------------------------------------------------------------------------------------------------------------------------

score(sin(x), sin(x), labels) # noiseless

score(y, y.corr, labels)
plot_series(x, y, y.corr)

score(y, y.acorr, labels)
plot_series(x, y, y.acorr)

score(y, y.uncorr, labels)
plot_series(x, y, y.acorr)


