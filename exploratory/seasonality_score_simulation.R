library(ggplot2)
library(lomb)
library(gridExtra)

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
    levels(labels) <- c("series 1", "series 2")
  }
  df <- data.frame(day = x, abundance = c(y1, y2), label = labels)
  p <- ggplot(df, aes(x = day, y = abundance, color = label)) +
    geom_point() +
    ylab("log abundance")
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
      0
    }
  })
  binary.y2 <- sapply(y2, function(yy) {
    if(yy >= mean.y2) {
      1
    } else {
      0
    }
  })
  n <- length(labels)
  season_v1 <- numeric(n)
  season_v1[labels == "wet"] <- 0
  season_v1[labels == "dry"] <- 1
  season_v2 <- numeric(n)
  season_v2[labels == "wet"] <- 1
  season_v2[labels == "dry"] <- 0
  match_score <- sum(binary.y1 == binary.y2)/n
  max_seasonality.y1 <- max(c(sum(binary.y1 == season_v1)/n,
                              sum(binary.y1 == season_v2)/n))
  max_seasonality.y2 <- max(c(sum(binary.y2 == season_v1)/n,
                              sum(binary.y2 == season_v2)/n))
  c(match_score, max_seasonality.y1, max_seasonality.y2)
}

## -----------------------------------------------------------------------------------------------------------------------------
##   Correlated series y, y.corr, anti-correlated series y.acorr, and uncorrelated series y.uncorr
## -----------------------------------------------------------------------------------------------------------------------------

noise_factor <- 0.1
x <- seq(from = 0, to = 10*pi, length.out = 1000)
y <- sin(x) + rnorm(1000)*noise_factor
y.corr <- sin(x) + rnorm(1000)*noise_factor
y.acorr <- sin(-x) + rnorm(1000)*noise_factor
y.uncorr1 <- rnorm(1000)*0.5
y.uncorr2 <- rnorm(1000)*0.5
labels <- as.factor(sin(x) > 0)
levels(labels) <- c("wet", "dry")
p1 <- plot_series(x, y, y.corr)
p2 <- plot_series(x, y, y.acorr)
p3 <- plot_series(x, y, y.uncorr1)
p4 <- plot_series(x, y.uncorr1, y.uncorr2)
grid.arrange(grobs = list(p1, p2, p3, p4), ncol = 2)

## -----------------------------------------------------------------------------------------------------------------------------
##   "Score" seasonal co-occurrence
## -----------------------------------------------------------------------------------------------------------------------------

plot_series(x, sin(x), sin(x))
score(sin(x), sin(x), labels) # noiseless

plot_series(x, y, y)
score(y, y, labels)

plot_series(x, y, y.corr)
score(y, y.corr, labels)

plot_series(x, y, y.acorr)
score(y, y.acorr, labels)

plot_series(x, y, y.uncorr1)
score(y, y.uncorr1, labels)

# zero-correlation (a coin flip) shows up as ~0.5 score
# anti-correlation is about 0
# positive correlation is about 1
plot_series(x, y.uncorr1, y.uncorr2)
score(y.uncorr1, y.uncorr2, labels)

