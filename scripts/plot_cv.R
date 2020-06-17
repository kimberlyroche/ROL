library(ggplot2)

cv_results <- read.table("output/cv_results.txt", header = FALSE, stringsAsFactors = TRUE, sep = "\t")
colnames(cv_results) <- c("value", "error_type", "host", "run_label")

# stick with log RMSE for now
error_type <- "log rmse"
cv_results <- cv_results[cv_results$error_type == error_type,]

# visualize error by model type ("run_label")
p <- ggplot(cv_results, aes(x = run_label, y = value)) +
  geom_boxplot() +
  xlab("Model type") +
  ylab("RMSE of log counts")
ggsave("runlabel_boxplot.png", p, units = "in", dpi = 150, height = 6, width = 8)

# visualize error by host for each model type
for(rl in unique(cv_results$run_label)) {
    p <- ggplot(cv_results[cv_results$run_label == rl,], aes(x = host, y = value)) +
                geom_boxplot() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                xlab("Host") +
                ylab(paste0("RMSE of log counts (model \"",rl,"\")"))
    ggsave(paste0("host_boxplot_",rl,".png"), p, units = "in", dpi = 150, height = 8, width = 12)
}

