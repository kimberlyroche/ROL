library(ggplot2)

label_run <- function(run_combo) {
  paste0("days",run_combo$days_decay,"_",
         "cov",run_combo$covariates,"_",
         "Gamma",run_combo$Gamma_scale,"_",
         "Xi",run_combo$Xi_scale)
}

cv_results <- read.table("output/cv_results.txt", header = FALSE, stringsAsFactors = TRUE, sep = "\t")
colnames(cv_results) <- c("value", "error_type", "host", "days_decay", "covariates", "Gamma_scale", "Xi_scale")
head(cv_results)

# stick with log RMSE for now
error_type <- "log rmse"
cv_results <- cv_results[cv_results$error_type == error_type,]

summ_df <- data.frame(error = c(), model = c())
run_combos <- unique(cv_results[,c("days_decay","covariates","Gamma_scale","Xi_scale")])
for(cidx in 1:nrow(run_combos)) {
  rc <- run_combos[cidx,]
  result_subset <- cv_results[(cv_results$days_decay == rc$days_decay &
                               cv_results$covariates == rc$covariates &
                               cv_results$Gamma_scale == rc$Gamma_scale &
                               cv_results$Xi_scale == rc$Xi_scale),]
  summ_df <- rbind(summ_df, data.frame(error = result_subset[result_subset$error_type == "log rmse",]$value,
                                       model = label_run(result_subset[1,])))
}

p <- ggplot(summ_df, aes(x = model, y = error)) +
  geom_boxplot() +
  xlab("Model type") +
  ylab("RMSE of log counts") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("runlabel_boxplot.png", p, units = "in", dpi = 150, height = 6, width = 8)

# visualize error by host for each model type
# p <- ggplot(cv_results[cv_results$run_label == rl,], aes(x = host, y = value)) +
#             geom_boxplot() +
#             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#             xlab("Host") +
#             ylab(paste0("RMSE of log counts (model \"",rl,"\")"))
# ggsave(paste0("host_boxplot_",rl,".png"), p, units = "in", dpi = 150, height = 8, width = 12)


