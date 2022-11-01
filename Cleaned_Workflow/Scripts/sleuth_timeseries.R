suppressMessages({
 library("sleuth")
})

suppressMessages({
 library("splines")
})

suppressMessages({
 library("dplyr")
})

#set input and output dirs
resultdir <- "/scratch/srb67793/2022VenusFlyTrap/sleuth"

setwd(resultdir)
#
# # ##################################
# # # No Prey Time Series
# # ##################################
#
# s2c <- read.csv("/scratch/srb67793/2022VenusFlyTrap/kallisto/Sleuth_NoPreyTimeSeries_SampleTissueConditionMinPath.csv", header = TRUE, stringsAsFactors = FALSE)
# s2c[] <- lapply(s2c, as.character)
# s2c$minutes <- as.numeric(s2c$minutes)
#
# time <- s2c$minutes
# full_design <- model.matrix(formula(~ ns(time, df = 4)))
# full_design
#
# so <- sleuth_prep(s2c, full_model = full_design)
# so <- sleuth_fit(so)
# so <- sleuth_fit(so, formula = ~ 1, fit_name = "reduced")
# so <- sleuth_lrt(so, "reduced", "full")
#
# pdf(file="SleuthNoPreyTimeSeriesQQPlot.pdf")
# plot_qq(so, test = 'reduced:full', test_type = 'lrt', sig_level = 0.05)
# dev.off()
#
# lrt_results <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
# sleuth_significant <- dplyr::filter(lrt_results, qval <= 0.05)
# write.csv(x = sleuth_significant, file = "SleuthNoPreyTimeSeriesResults_qval_0.05.csv", row.names = FALSE)
# write.csv(x = lrt_results, file = "SleuthNoPreyTimeSeriesResults.csv", row.names = FALSE)
# pdf(file="SleuthNoPreyTimeSeriesTop20HeatMap.pdf")
# plot_transcript_heatmap(so, head(lrt_results, n = 20)$target_id, 'est_counts')
# dev.off()
#
# pdf(file="SleuthNoPreyTimeSeriesResults.pdf")
# for(i in sleuth_significant$target_id[1:151]) {
#   p1 <- plot_bootstrap(sleuth_object, i, units = "tpm", color_by = "condition")
#   print(p1)
# }
# dev.off()

##################################
# Prey Time Series
##################################

s2c <- read.csv("/scratch/srb67793/2022VenusFlyTrap/kallisto/Sleuth_PreyTimeSeries_SampleTissueConditionMinPath.csv", header = TRUE, stringsAsFactors = FALSE)
s2c[] <- lapply(s2c, as.character)
s2c$minutes <- as.numeric(s2c$minutes)

time <- s2c$minutes
full_design <- model.matrix(formula(~ ns(time, df = 4)))
full_design

so <- sleuth_prep(s2c, full_model = full_design)
so <- sleuth_fit(so)
so <- sleuth_fit(so, formula = ~ 1, fit_name = "reduced")
so <- sleuth_lrt(so, "reduced", "full")

pdf(file="SleuthPreyTimeSeriesQQPlot.pdf")
plot_qq(so, test = 'reduced:full', test_type = 'lrt', sig_level = 0.05)
dev.off()

lrt_results <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
sleuth_significant <- dplyr::filter(lrt_results, qval <= 0.05)
write.csv(x = sleuth_significant, file = "SleuthPreyTimeSeriesResults_qval_0.05.csv", row.names = FALSE)

pdf(file="SleuthPreyTimeSeriesTop20HeatMap.pdf")
plot_transcript_heatmap(so, head(lrt_results, n = 20)$target_id, 'est_counts')
dev.off()

pdf(file="SleuthNoPreyTimeSeriesResults.pdf")
for(i in sleuth_significant$target_id[1:880]) {
  p1 <- plot_bootstrap(so, i, units = "tpm", color_by = "condition")
  print(p1)
}
dev.off()

quit(save="no")

#
# tmp <- so$obs_raw %>% dplyr::filter(target_id == 'INPUT TARGET ID')
# tmp <- dplyr::full_join(so$sample_to_covariates, tmp, by = 'sample')
# tmp
# tmp <- transform(tmp, time = as.numeric(time))
# ggplot(tmp, aes(x=time, y=est_counts))
#   + geom_point(shape=1)
#   + geom_smooth(method = loess)
#
# pdf(file="SleuthTimeSeriesTop20HeatMap.pdf")
# print(tmp)
# dev.off()
