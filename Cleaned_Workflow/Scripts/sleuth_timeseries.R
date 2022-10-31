suppressMessages({
 library("sleuth")
})

suppressMessages({
 library("splines")
})

#set input and output dirs
datapath <- "/scratch/srb67793/2022VenusFlyTrap/kallisto/quant/"
resultdir <- "/scratch/srb67793/2022VenusFlyTrap/sleuth"

setwd(resultdir)

s2c <- read.table("VFT_samples_condition_time_path.csv", header = TRUE, stringsAsFactors = FALSE)
s2c <- dplyr::mutate(s2c, path = paths)
s2c[] <- lapply(s2c, as.character)
s2c

day <- s2c$time
full_design <- model.matrix(formula(~ ns(time, df = 4)))
full_design

so <- sleuth_prep(s2c, full_model = full_design, target_mapping = t2g)
so <- sleuth_fit(so)
so <- sleuth_fit(so, formula = ~ 1, fit_name = "reduced")
so <- sleuth_lrt(so, "reduced", "full")

pdf(file="SleuthTimeSeriesQQPlot.pdf")
plot_qq(so, test = 'reduced:full', test_type = 'lrt', sig_level = 0.05)
dev.off()

lrt_results <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
table(lrt_results[,"qval"] < 0.05)

lrt_results %>% head(n = 20) %>% dplyr::select(target_id, qval)
write.csv(x = lrt_results, file = "SleuthTimeSeriesResults.csv", row.names = FALSE)

pdf(file="SleuthTimeSeriesTop20HeatMap.pdf")
plot_transcript_heatmap(so, head(lrt_results, n = 20)$target_id, 'est_counts')
dev.off()
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
