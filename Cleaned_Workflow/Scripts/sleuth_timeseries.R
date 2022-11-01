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

# read data into sleuth_object, import bootstrap summaries and TPMs, and perform normalization/filtering steps
sleuth_object <- sleuth_prep(s2c, full_model = full_design, extra_bootstrap_summary = TRUE, read_bootstrap_tpm=TRUE)

# estimate parameters for the full linear model that includes the conditions as factors
sleuth_object <- sleuth_fit(sleuth_object)

# estimate parameters for the reduced linear model that assumes equal transcript abundances in both conditions
sleuth_object <- sleuth_fit(sleuth_object, ~1, 'reduced')

# perform likelihood ratio test to identify transcripts whose fit is significantly better under full model relative to reduced model
sleuth_object <- sleuth_lrt(sleuth_object, 'reduced', 'full')

# check that sleuth object is OK
models(sleuth_object)

#summarize the sleuth results for DE genes with q-values < 0.05
sleuth_table <- sleuth_results(sleuth_object, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
write.csv(x = sleuth_significant, file = "SleuthPreyTimeSeriesResults_qval_0.05.csv", row.names = FALSE)

pdf(file="SleuthPreyTimeSeriesQQPlot.pdf")
plot_qq(sleuth_object, test = 'reduced:full', test_type = 'lrt', sig_level = 0.05)
dev.off()

pdf(file="SleuthPreyTimeSeriesTop20HeatMap.pdf")
plot_transcript_heatmap(sleuth_object, head(sleuth_significant, n = 20)$target_id, 'est_counts')
dev.off()

pdf(file="SleuthPreyTimeSeriesResults.pdf")
for(i in sleuth_significant$target_id[1:881]) {
  p1 <- plot_bootstrap(sleuth_object, i, units = "tpm")
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
